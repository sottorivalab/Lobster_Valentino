#!/bin/bash

# use R version 4.2.3

# post-processing code
user_dir="/scratch/madeleine.dale/lobster_paper_fresh/lobster4publishing" # replace with path for results directory

#Rscript install_required_R_packages.R

# create required subdirectories
#cd ${user_dir}
#Rscript vignettes/post_processing_VALE/create_required_subdirectories.R
#cd -

Rscript denovo_clustering_VALE.R 
Rscript annotate_ploidy_mix_VALE.R 
Rscript scatterplots_VALE.R 
Rscript merge_VALE_dataframes.R 

# unite scatterplot pdfs 
cd ${user_dir}/results/vale/vale_scatterplots_pass/clinical_scatters/significance
pdfunite $(cat ${user_dir}/vale/scatterplot_file_list.txt) merged.pdf

cd ${user_dir}/results/vale/vale_scatterplots_pass/clinical_scatters/significance_untransformed
pdfunite $(cat ${user_dir}/vale/scatterplot_file_list.txt) merged.pdf
