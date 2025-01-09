#!/bin/bash

# Check if input and output arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_fasta> <output_fasta>"
    exit 1
fi

# Input and output file paths
input="$1"
output="$2"

seq=""
header=""

{
    while IFS= read -r line; do
        if [[ ${line:0:1} == ">" ]]; then
            if [[ -n $seq && ${#seq} -ge 120 ]]; then
                echo "$header"
                echo "$seq"
            fi
            header="$line"
            seq=""  # Reset sequence
        else
            seq+="$line"  # Concatenate sequence lines
        fi
    done < "$input" > "$output"

    # Check the last sequence
    if [[ -n $seq && ${#seq} -ge 120 ]]; then
        echo "$header"
        echo "$seq"
    fi
}
