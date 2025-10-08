library(ggplot2)
library(ggrepel)
library(plyr)
library(dplyr)
library(R.utils)
library(mobster)
library(lubridate)
library(tidyr)
library(stats)
library(DEoptim)
library(abind)
library(combinat)
library(vcfR)
library(purrr)

##########################################
run_lobster <- function(df_not_used,patient_str,cohort,filtering_str,user_dir){
  
  xstar = 0.05
  Kvec = seq(3) # max no. of possible beta components
  
  for (ploidy in names(df_not_used)){
    if (ploidy!='ploidy_mix'){
        
      df = readRDS(sprintf('%s/results/vale/vale_fits_pass_original/transformed_%s_%s.Rds',user_dir,patient_str,ploidy))
      
      if (!is.null(df)){
        df = df$best$data
        # this clustering info will be recalculated by lobster
        df$cluster = NULL
        df$kmeans_cluster = NULL
         
        mfit = df %>% mobster_fit(K=Kvec, samples=1, parallel = TRUE,N_cutoff=1,maxIter = 15, genSA_steps = 3, quantile=TRUE,log_MMscaled = FALSE, untf_convolved = TRUE, xstar = xstar, Gopt = 'DE', no_filtering=TRUE)
        
        saveRDS(mfit, sprintf('%s/results/%s/%s_fits_from_dataframe/transformed_%s_%s.Rds',user_dir,cohort,cohort,patient_str,ploidy))
  }}}
}
################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Provide a patient no. as the argument to this program.", call. = FALSE)
}
i <- as.integer(args[1])
print(i)

config = mobster:::return_config_dict()
user_dir = config$paths$user_dir

cohort = 'vale' # lower case
filtering='pass'
filtering_str='pass'

df = mobster:::format_vale_wes(i, filtering)
patient_str = df$patient_split[[1]]$patient[[1]] 
df = mobster:::split_ploidies(df)
run_lobster(df,patient_str,cohort,filtering_str,user_dir)
      