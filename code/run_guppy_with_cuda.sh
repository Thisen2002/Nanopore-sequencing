#!/bin/bash

# Check if correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

dir_input=$1
dir_output=$2

# Create the output/raw_data directory if it doesn't exist
mkdir -p "$dir_output"

# Flow cell and kit can be overridden to match the run, e.g. KIT=SQK-LSK110
flowcell="${FLOWCELL:-FLO-MIN106}"
kit="${KIT:-SQK-PCB109}"

guppy_basecaller -i "$dir_input" -s "$dir_output" --flowcell "$flowcell" --kit "$kit" -x 'cuda:all'

