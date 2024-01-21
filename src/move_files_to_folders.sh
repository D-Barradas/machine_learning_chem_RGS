#!/bin/bash

# Set the directory to process
dir_to_process="../data"

# Initialize the counter and folder number
counter=0
folder_num=1

# Loop through all files in the directory
for file in `cat $dir_to_process/compunds_cids.txt`; do

  # Check if the file ends with ".csv"
  if [ -f $dir_to_process/preliminar_data_${file}.csv ]  ; then
    # If the counter reaches 1000, create a new folder and reset the counter
    if [ $counter -eq 1000 ]; then
      mkdir -p "$dir_to_process/folder_$folder_num"
      counter=0
      folder_num=$((folder_num + 1))
    fi

    # Move the file to the current folder
    mv "$dir_to_process/preliminar_data_${file}.csv" "$dir_to_process/folder_$folder_num"/

    # Increment the counter
    counter=$((counter + 1))
  fi
done

