#!/bin/bash

# Change to the build directory
cd "$(dirname "$0")/build"

# Define the folders to process
folders=("Cordeau" "Solomon")

# Define the values for m and s
m_values=(1 2 3 4)
s_values=(1 2 3 4)

# Loop over the folders
for folder in "${folders[@]}"
do
  # Find all files in the folder
  files=$(find "./instances/$folder" -name '*.txt' -type f)
  
  # Loop over the files
  for file in $files
  do
    # Remove the path prefix from the file name
    file=$(basename "$file")
    
    # Loop over the values of m and s
    for m in "${m_values[@]}"
    do
      for s in "${s_values[@]}"
      do
        # Construct the command to run the binary with the current arguments
        cmd="./AMTOPTW -f $folder -i ${file%.*} -m $m -s $s"
        
        # Print the command being executed
        echo "Running command: $cmd"
        
        # Run the binary with the current arguments
        $cmd
      done
    done
  done
done
