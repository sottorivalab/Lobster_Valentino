#!/bin/bash
#SBATCH --job-name=rscripts
#SBATCH --partition=cpuq
#SBATCH --ntasks=3       
#SBATCH --output=rscript_%j.out
#SBATCH --mem-per-cpu=50000M
#SBATCH --time=3-0:00:00

module load R/4.2.3
module load R/rstudio-dependencies
module load R/4.2.3-rstudio-dependencies
Rscript --verbose run_lobster_i_dataframe.R ${i} &> patient_vale_dataframe_"$i"_job"$SLURM_JOB_ID".out

