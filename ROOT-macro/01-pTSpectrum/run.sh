#!/bin/bash

macro=$1

# Check if the file exists
if [[ ! -f "$macro" ]]; then
    echo "Error: $macro not found!"
    exit 1
fi

# Run the ROOT macro in batch mode
root -l <<EOF
.L $1
read_particles()
EOF
