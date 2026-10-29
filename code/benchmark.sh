#!/bin/bash

# Run the pipeline stage by stage and record runtime, CPU, memory, GPU and read statistics for each stage

# Check if correct number of arguments are provided
if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage: $0 <input_directory> <output_directory> [--fastq-input]"
    echo "  --fastq-input  skip basecalling; <input_directory> holds .fastq/.fastq.gz reads"
    echo "  Set FLOWCELL and KIT to match the run, e.g. KIT=SQK-LSK110 $0 data/ out/"
    exit 1
fi

dir_input=$1
dir_output=$2
fastq_input=false
if [ "$#" -eq 3 ]; then
    if [ "$3" == "--fastq-input" ]; then
        fastq_input=true
    else
        echo "Invalid argument: $3"
        exit 1
    fi
fi

script_dir=$(cd "$(dirname "$0")" && pwd)
dir_raw_data="$dir_output/raw_data"
dir_fastq_data="$dir_output/uncleaned_fastq_output"
dir_cleaned_fastq_data="$dir_output/cleaned_fastq_output"
dir_fasta_data="$dir_output/fasta_output"
dir_logs="$dir_output/benchmark_logs"
report="$dir_output/benchmark_report.tsv"
summary="$dir_output/benchmark_summary.txt"

# Check required tools before starting
required="porechop seqtk zip python3 /usr/bin/time"
if [ "$fastq_input" = false ]; then
    required="guppy_basecaller $required"
fi
missing=""
for tool in $required; do
    command -v "$tool" > /dev/null || missing="$missing $tool"
done
if [ -n "$missing" ]; then
    echo "Missing required tools:$missing"
    exit 1
fi

# Clear previous benchmark outputs only (not the whole output directory)
rm -rf "$dir_raw_data" "$dir_fastq_data" "$dir_cleaned_fastq_data" "$dir_fasta_data" "$dir_logs"
mkdir -p "$dir_logs"
echo -e "stage\texit_code\twall_s\tcpu_s\tcpu_pct\tmax_rss_mb\tgpu_util_avg_pct\tgpu_mem_max_mb" > "$report"

# Run one stage under /usr/bin/time, sampling GPU usage while it runs
run_stage() {
    local name=$1
    shift
    local gpu_log="$dir_logs/$name.gpu.csv"
    local gpu_pid=""

    echo
    echo "=== Stage: $name"
    if command -v nvidia-smi > /dev/null; then
        nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv,noheader,nounits -l 1 > "$gpu_log" 2>/dev/null &
        gpu_pid=$!
    fi

    /usr/bin/time -v -o "$dir_logs/$name.time" "$@" > "$dir_logs/$name.log" 2>&1
    local rc=$?

    if [ -n "$gpu_pid" ]; then
        kill "$gpu_pid" 2>/dev/null
        wait "$gpu_pid" 2>/dev/null
    fi

    local wall cpu pct rss gpu_util gpu_mem
    wall=$(awk -F': ' '/Elapsed \(wall clock\)/ {n=split($2,t,":"); s=0; for(i=1;i<=n;i++) s=s*60+t[i]; printf "%.2f", s}' "$dir_logs/$name.time")
    cpu=$(awk -F': ' '/User time|System time/ {s+=$2} END {printf "%.2f", s}' "$dir_logs/$name.time")
    pct=$(awk -F': ' '/Percent of CPU/ {sub("%","",$2); print $2}' "$dir_logs/$name.time")
    rss=$(awk -F': ' '/Maximum resident set size/ {printf "%.0f", $2/1024}' "$dir_logs/$name.time")
    gpu_util="NA"
    gpu_mem="NA"
    if [ -s "$gpu_log" ]; then
        gpu_util=$(awk -F', ' '{s+=$1; n++} END {printf "%.1f", s/n}' "$gpu_log")
        gpu_mem=$(awk -F', ' '$2>m {m=$2} END {print m+0}' "$gpu_log")
    fi

    echo -e "$name\t$rc\t$wall\t$cpu\t$pct\t$rss\t$gpu_util\t$gpu_mem" >> "$report"
    echo "exit=$rc wall=${wall}s cpu=${cpu}s max_rss=${rss}MB gpu_util=${gpu_util}% (log: $dir_logs/$name.log)"
    if [ "$rc" -ne 0 ]; then
        echo "Stage $name failed, see $dir_logs/$name.log"
        exit "$rc"
    fi
}

# Copy .fast5 files the same way main.sh does (a script string, since /usr/bin/time cannot run shell functions)
collect_fast5='
    mkdir -p "$2"
    find "$1" \( -name "pass_barcode*" -o -name "fail_barcode*" -o -name "fast5_pass*" -o -name "fast5_fail*" \) -type d \
        -exec sh -c "for dir do cp -t \"\$0\" -- \"\$dir\"/*.fast5; done" "$2" {} +
    ls "$2"/*.fast5 > /dev/null 2>&1 || { echo "No .fast5 files found in $1"; exit 1; }
'

if [ "$fastq_input" = true ]; then
    pass_reads="$dir_input"
    fail_reads=""
else
    run_stage 1_collect_fast5 bash -c "$collect_fast5" collect "$dir_input" "$dir_raw_data"
    run_stage 2_basecalling "$script_dir/run_guppy_with_cuda.sh" "$dir_raw_data" "$dir_fastq_data"
    pass_reads="$dir_fastq_data/pass"
    fail_reads="$dir_fastq_data/fail"
fi
run_stage 3_porechop "$script_dir/run_data_filter.sh" "$pass_reads" "$dir_cleaned_fastq_data"
run_stage 4_fastq_to_fasta "$script_dir/fastq_to_fasta.sh" "$dir_cleaned_fastq_data" "$dir_fasta_data"

# Read statistics at each point of the pipeline
echo
echo "=== Computing read statistics"
stats_header="dataset\tfiles\treads\tbases\tmean_len\tn50\tmax_len\tmean_q"
stats_table=$(
    echo -e "$stats_header"
    echo -e "basecalled_pass\t$(python3 "$script_dir/read_stats.py" "$pass_reads")"
    if [ -n "$fail_reads" ]; then
        echo -e "basecalled_fail\t$(python3 "$script_dir/read_stats.py" "$fail_reads")"
    fi
    echo -e "porechop_cleaned\t$(python3 "$script_dir/read_stats.py" "$dir_cleaned_fastq_data")"
    echo -e "final_fasta\t$(python3 "$script_dir/read_stats.py" "$dir_fasta_data/converted_data.fasta")"
)

# Derived scores
scores=$(echo "$stats_table" | awk -F'\t' '
    {reads[$1]=$3; bases[$1]=$4}
    END {
        if ("basecalled_fail" in reads) {
            total = reads["basecalled_pass"] + reads["basecalled_fail"]
            if (total > 0) printf "Basecalling pass rate:        %.1f%% of reads\n", 100*reads["basecalled_pass"]/total
        }
        if (reads["basecalled_pass"] > 0) {
            printf "Porechop read retention:      %.1f%%\n", 100*reads["porechop_cleaned"]/reads["basecalled_pass"]
            printf "Porechop base retention:      %.1f%%\n", 100*bases["porechop_cleaned"]/bases["basecalled_pass"]
        }
    }')
guppy_rate=$(grep -ho 'samples/s: [0-9.e+]*' "$dir_fastq_data"/guppy_basecaller_log-*.log 2>/dev/null | tail -1)
total_wall=$(awk -F'\t' 'NR>1 {s+=$3} END {printf "%.1f", s}' "$report")
basecall_wall=$(awk -F'\t' '$1=="2_basecalling" {print $3}' "$report")

{
    echo "Pipeline benchmark - $(date)"
    echo "Host: $(hostname)  CPUs: $(nproc)  RAM: $(free -h | awk '/Mem:/ {print $2}')"
    echo "CPU:  $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //')"
    if command -v nvidia-smi > /dev/null; then
        echo "GPU:  $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader | paste -sd ';')"
    fi
    echo "Input: $dir_input  FLOWCELL=${FLOWCELL:-FLO-MIN106} KIT=${KIT:-SQK-PCB109}"
    echo "Tools: $(guppy_basecaller --version 2>/dev/null | grep -m1 -o 'Version [^ ]*') porechop $(porechop --version 2>/dev/null) seqtk $(seqtk 2>&1 | awk '/Version/ {print $2}')"
    echo
    echo "--- Runtime and resources per stage"
    column -t -s $'\t' "$report"
    echo
    echo "Total wall time: ${total_wall}s"
    [ -n "$guppy_rate" ] && echo "Guppy throughput: $guppy_rate"
    if [ -n "$basecall_wall" ]; then
        awk -v r="$(echo "$stats_table" | awk -F'\t' '$1 ~ /^basecalled_/ {s+=$3} END {print s}')" -v w="$basecall_wall" \
            'BEGIN {if (w > 0) printf "Basecalling speed: %.1f reads/s\n", r/w}'
    fi
    echo
    echo "--- Read statistics"
    echo "$stats_table" | column -t -s $'\t'
    echo
    echo "--- Scores"
    echo "$scores"
} > "$summary"

echo
cat "$summary"
echo
echo "Saved: $report and $summary"
