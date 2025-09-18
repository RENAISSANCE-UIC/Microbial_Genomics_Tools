#!/bin/bash

# Set backup directory
BACKUP_DIR="$HOME/conda_env_backups/$(date +%Y%m%d)"
mkdir -p "$BACKUP_DIR"

echo "Backing up conda environments to $BACKUP_DIR..."

# Get list of conda environments (excluding base if desired)
ENVS=$(conda env list | grep -v "^#" | awk '{print $1}')

# Loop through each environment and export it
for ENV in $ENVS; do
    # Skip empty lines
    if [ -z "$ENV" ]; then
        continue
    fi
    
    echo "Exporting environment: $ENV"
    
    # Export the environment to YAML
    conda env export -n "$ENV" > "$BACKUP_DIR/${ENV}_environment.yml"
    
    # Alternative: export only packages directly installed by user
    # conda env export -n "$ENV" --from-history > "$BACKUP_DIR/${ENV}_environment_minimal.yml"
    
    echo "Completed export of $ENV"
done
