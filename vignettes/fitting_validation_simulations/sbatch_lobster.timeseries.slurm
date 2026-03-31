#!/bin/bash
#SBATCH --job-name=rscripts
#SBATCH --partition=cpuq
#SBATCH --ntasks=3       
#SBATCH --output=logs/rscript_%j.out
#SBATCH --mem-per-cpu=5000M
#SBATCH --time=24:00:00

module load R/4.2.3
module load R/rstudio-dependencies
module load R/4.2.3-rstudio-dependencies

FIT_STR=$1
Rscript --verbose run_lobster_timeseries_simulation.R "$FIT_STR" &> logs/lobster_validation_test_job"$SLURM_JOB_ID".out

