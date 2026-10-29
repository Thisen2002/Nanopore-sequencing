#!/usr/bin/env python3
"""Print read statistics for FASTQ/FASTA files (plain or gzipped) as one tab-separated line.

Usage: read_stats.py <file_or_directory> [...]
Columns: files, reads, bases, mean_len, n50, max_len, mean_q
mean_q is the average of per-read mean Q scores (computed from error
probabilities, as Guppy does); it is NA for FASTA input.
"""
import gzip
import math
import sys
from pathlib import Path

SUFFIXES = (".fastq", ".fq", ".fasta", ".fa", ".fastq.gz", ".fq.gz", ".fasta.gz", ".fa.gz")
ERROR_PROB = [10 ** (-q / 10) for q in range(94)]


def open_reads(path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def read_records(path):
    """Yield (length, mean_q or None) for every read in the file."""
    with open_reads(path) as handle:
        first = handle.read(1)
    if not first:
        return
    with open_reads(path) as handle:
        if first == "@":
            while True:
                header = handle.readline()
                if not header:
                    break
                seq = handle.readline().rstrip("\n")
                handle.readline()
                qual = handle.readline().rstrip("\n")
                if not seq:
                    yield 0, None
                    continue
                mean_err = sum(ERROR_PROB[ord(c) - 33] for c in qual) / len(qual)
                yield len(seq), -10 * math.log10(mean_err)
        else:
            length = None
            for line in handle:
                if line.startswith(">"):
                    if length is not None:
                        yield length, None
                    length = 0
                else:
                    length += len(line.strip())
            if length is not None:
                yield length, None


def collect_files(args):
    files = []
    for arg in args:
        path = Path(arg)
        if path.is_dir():
            files.extend(sorted(p for p in path.iterdir() if p.name.endswith(SUFFIXES)))
        elif path.is_file():
            files.append(path)
    return files


def main():
    files = collect_files(sys.argv[1:])
    lengths, qualities = [], []
    for path in files:
        for length, q in read_records(path):
            lengths.append(length)
            if q is not None:
                qualities.append(q)

    if not lengths:
        print(f"{len(files)}\t0\t0\tNA\tNA\tNA\tNA")
        return

    bases = sum(lengths)
    running, n50 = 0, 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= bases / 2:
            n50 = length
            break
    mean_q = f"{sum(qualities) / len(qualities):.2f}" if qualities else "NA"
    print(f"{len(files)}\t{len(lengths)}\t{bases}\t{bases / len(lengths):.0f}\t{n50}\t{max(lengths)}\t{mean_q}")


if __name__ == "__main__":
    main()
