In order to run the Lobster tool and post-processing for the VALE cohort:

The relevant scripts are in the vignettes/processing_VALE directory.

Steps
-----
- Build the package as a project in Rstudio (R version 4.2.3)

- cd to vignettes/processing_VALE directory

- run the script "./slurm_run_command.sh" for the Lobster fitting
  |->  this creates a slurm job for each patient
  |->  resulting fit objects are placed inside {user_dir}/results/vale/vale_fits_pass
  
- run the script "./run_postprocessing.sh" to run the post-processing
  |->  this will output merged pdfs of the "transformed" and "untransformed" scatterplots
       (paths {user_dir}/results/vale/vale_scatterplots_pass/clinical_scatters/{significance|significance_transformed}/merged.pdf)
       as well as a dataframe containing all the Lobster-annotated mutations over the cohort 
       ({user_dir}/results/vale/vale_pass_ploidy_mix_new/merged_final_BH_0.05.Rds)
       
These scripts contain a "user_dir" results directory path that needs to be provided by the user;
the scripts should then create the required directory structure

The user also needs to put the paths to the maf and ichor seg directories and to the user_dir in the user_config.yml file in the project home directory.
These maf and ichor seg directories can be downloaded from the Mendeley page. 

The fits used in the paper are given in results/vale_fits_pass_original.

Please address any questions or problems in running the code to madeleine.dale@fht.org.

------------
Copyright Madeleine Dale October 2025

This code is adapted from the MOBSTER package published in: 
Caravagna G, Sanguinetti G, Graham TA, Sottoriva A. The MOBSTER R package for tumour subclonal deconvolution from bulk DNA whole-genome sequencing data. BMC Bioinformatics. Springer Science and Business Media LLC; 2020;21:531.

It uses the DEoptim global optimisation package published in:
Mullen K, Ardia D, Gil D, Windover D, Cline J. DEoptim: AnRPackage for global optimization by differential evolution. J Stat Softw. Foundation for Open Access Statistic; 2011;40:1–26.

