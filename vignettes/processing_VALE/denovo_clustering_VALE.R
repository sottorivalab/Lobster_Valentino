library(cowplot)
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
library(stringr)
library(yaml)

config = mobster:::return_config_dict()
user_dir = config$paths$user_dir
from_dataframe = config$from_dataframe

cohort = 'VALE'
vale_df = read.csv(sprintf('%s/vale/VALENTINO_samplesheet.csv', user_dir)) 
imax = 37
filtering = 'pass'
filtering_str='pass'

results_dir = sprintf('%s/results',user_dir)
if (from_dataframe==TRUE){
  fit_path = sprintf('%s/%s/%s_fits_from_dataframe',results_dir,tolower(cohort),tolower(cohort))
} else {
  fit_path = sprintf('%s/%s/%s_fits_%s',results_dir,tolower(cohort),tolower(cohort),filtering_str)
}

ploidies = c('NA',as.character(seq(5)))

sample_names = vale_df$sample
patient_names = vale_df$case

for (i in seq(1,imax)){
  patient = unique(patient_names)[[i]]
  full_df = mobster:::format_vale_wes(i, filtering)
  df = mobster:::split_ploidies(full_df)
  full_df = mobster:::append_VAF_time(full_df$patient_time_split,1) 
  
  # create new ploidy_mix
  saveRDS(df[['ploidy_mix']], sprintf('%s/%s/%s_%s_ploidy_mix_new/%s_%s_ploidy_mix.Rds',results_dir,tolower(cohort),tolower(cohort),filtering_str,patient,filtering_str))
  
  full_df$VAF = full_df$VAF %>% mutate_all(~ ifelse(is.nan(.), NA, .))
  print(patient)
  print(paste0('rows before excluding all NA: ',nrow(full_df)))
  full_df <- full_df[apply(full_df$VAF, 1, function(x) !all(is.na(unlist(x)))), ]
  print(paste0('rows after excluding all NA: ',nrow(full_df)))
  total_rows = nrow(full_df)
  
  for (ploidy_idx in seq(1,length(ploidies))){
    ploidy = ploidies[ploidy_idx]
  
    # generate un-sampled dataframe
    if (ploidy!='NA'){
      full_df_idx = which(apply(full_df$ploidy, 1, function(row) all(row == as.numeric(ploidy))))
      full_df_ploided = full_df[full_df_idx,]
      print(sprintf('nrow(full_df_ploided): %d',nrow(full_df_ploided)))
    } else { full_df_ploided = full_df } 
    
    try({
        # ploidies for which fit doesn't exist will be caught by annotate_ploidy_mix.Rmd
        mfit = readRDS(sprintf('%s/transformed_%s_ploidy_%s.Rds',fit_path,patient,ploidy))
        mfit[['mutations_over_all_ploidies']] = total_rows
        saveRDS(mfit,sprintf('%s/annotated/transformed_%s_ploidy_%s_%s_denovo.Rds',fit_path,unique(patient_names)[i],ploidy,filtering_str))
        
        selected_fit = mobster:::select_fit(mfit,patient,ploidy,cohort,user_dir)
        mfit$annotated = mobster:::Clusters_denovo(selected_fit, full_df_ploided)
        saveRDS(mfit,sprintf('%s/annotated/transformed_%s_ploidy_%s_%s_denovo.Rds',fit_path,unique(patient_names)[i],ploidy,filtering_str))
    })
  } 
}