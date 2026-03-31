#!/bin/bash

# Directory containing input Rds files
INPUT_DIR="/scratch/madeleine.dale/lobster_paper_fresh/lobster4publishing.VALE/vignettes/fitting_validation/sampling_simulations/dataframes"

if [ ! -d "$INPUT_DIR" ]; then
  echo "Error: Directory $INPUT_DIR not found."
  exit 1
fi

for file in "$INPUT_DIR"/*.Rds; do
    filename=$(basename "$file")
    id="${filename%.*}"

    echo "Submitting job for ID: $id"
    
    sbatch --job-name="lob_$id" --output="logs/lob_$id.out" --error="logs/lob_$id.err" sbatch_lobster.single_timepoint.slurm "$id"
done