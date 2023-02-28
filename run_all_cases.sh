#!/bin/bash

#Change to the build directory
cd "$(dirname "$0")/build"

#Define the folders to process
folders=("Cordeau" "Solomon")

#Define the values for m and s
m_values=(1 2 3 4)
s_values=(1 2 3 4)

#Create the header row for the excel file
echo "instance,m,s,score,time" > benchmarks.csv

#Loop over the values of m and s
for m in "${m_values[@]}"
do
    for s in "${s_values[@]}"
    do
        # Loop over the folders
        for folder in "${folders[@]}"
        do
            # Find all files in the folder and sort them by name
            files=$(find "./instances/$folder" -name '*.txt' -type f | sort)

            # Loop over the files
            for file in $files
            do
                # Remove the path prefix from the file name
                file=$(basename "$file")
                
                # Construct the command to run the binary with the current arguments
                cmd="./AMTOPTW -f $folder -i ${file%.*} -m $m -s $s"

                # Print the command being executed
                echo "Running command: $cmd"

                # Run the binary with the current arguments and capture the output
                output=$($cmd)

                # Extract the score and time from the output
                score=$(echo "$output" | grep "Score" | awk '{print $2}')
                time=$(echo "$output" | grep "Time" | awk '{print $2}')

                # Write the results to the excel file
                echo "${file%.*},$m,$s,$score,$time" >> benchmarks.csv
            done
        done
    done
done