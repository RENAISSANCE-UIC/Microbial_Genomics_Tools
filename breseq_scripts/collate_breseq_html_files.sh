
#!/bin/bash

# Define the base directory
BASE_DIR="/home/weackerm/Desktop/shortcut-to-shared-nas/Stabryla_Data/Motility_AMR/results/breseq"

# Change to the base directory
cd "$BASE_DIR" || { echo "Failed to change directory to $BASE_DIR"; exit 1; }

# Loop through each subdirectory in the base directory
for folder in "$BASE_DIR"/*; do
    if [ -d "$folder" ]; then
        folder_name=$(basename "$folder")
        source_file="$folder/output/summary.html"
        destination_file="$BASE_DIR/${folder_name}_summary.html"

        if [ -f "$source_file" ]; then
            cp "$source_file" "$destination_file"
            echo "Copied and renamed $source_file to $destination_file"
        else
            echo "File $source_file does not exist"
        fi
    fi
done


