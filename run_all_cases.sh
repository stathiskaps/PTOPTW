#!/bin/bash

# Change to the build directory
cd "$(dirname "$0")/build"

# Get folder name and file name from command line arguments
folder=$1
file=$2

# Define values of s and m
s_values=(1 2 3 4)
m_values=(1 2 3 4)

# Loop over all combinations of s and m
for s in "${s_values[@]}"
do
    for m in "${m_values[@]}"
    do
        # Call the binary with the given inputs
        ./AMTOPTW -f $folder -i $file -m $m -s $s -r
    done
done