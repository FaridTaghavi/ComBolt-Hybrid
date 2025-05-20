#!/bin/bash

# Check if a filename is given
if [ $# -ne 1 ]; then
    echo "Usage: $0 input_file"
    exit 1
fi

input_file="$1"
output_prefix="chunk"
chunk_num=0
outfile=""

while IFS= read -r line || [ -n "$line" ]; do
    # Count number of fields in the line
    num_fields=$(echo "$line" | awk '{print NF}')

    if [ "$num_fields" -eq 7 ]; then
        # Start a new chunk
        chunk_num=$((chunk_num + 1))
        outfile="${output_prefix}_${chunk_num}.txt"
        echo "$line" > "$outfile"
    else
        # Append to current chunk
        if [ -n "$outfile" ]; then
            echo "$line" >> "$outfile"
        fi
    fi
done < "$input_file"

echo "Split into $chunk_num chunks."

