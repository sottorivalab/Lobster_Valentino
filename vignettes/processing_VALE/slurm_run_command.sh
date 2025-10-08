#!/bin/bash

user_dir="/scratch/madeleine.dale/lobster_paper_fresh/lobster4publishing" # replace with path for results directory

Rscript ${user_dir}/vignettes/post_processing_VALE/install_required_R_packages.R

# create required subdirectories
cd ${user_dir}
Rscript ${user_dir}/vignettes/post_processing_VALE/create_required_subdirectories.R
cd -

for i in $(seq 1 37); do 
    sbatch --export=ALL,i=$i sbatch_script_vale.sh;
done
