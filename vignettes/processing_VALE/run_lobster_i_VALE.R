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
run_lobster <- function(df,patient_str,cohort,filtering_str,user_dir){
  
  xstar = 0.05
  Kvec = seq(3) # max no. of possible beta components
  
  for (ploidy in names(df)){
    if (ploidy!='ploidy_mix'){
  
      if (nrow(df[[ploidy]])>0){
        mfit = df[[ploidy]] %>% mobster_fit(K=Kvec, samples=1, parallel = TRUE,N_cutoff=1,maxIter = 15, genSA_steps = 3, quantile=TRUE,log_MMscaled = FALSE, untf_convolved = TRUE, xstar = xstar, Gopt = 'DE')
        
        saveRDS(mfit, sprintf('%s/results/%s/%s_fits_%s/transformed_%s_%s.Rds',user_dir,cohort,cohort,filtering_str,patient_str,ploidy))
         
        plot_list = plot(mfit$best,log_MMscaled = FALSE, untf = FALSE, transform_axis_back = TRUE, highlight_NA = FALSE, annotate_genes = FALSE)
        
        for (i in 1:length(plot_list)){
          ggsave(plot = plot_list[[i]], filename = sprintf('%s/results/%s/%s_plots_%s/transformed_%s_%s_t%d.pdf',user_dir,cohort,cohort,filtering_str,patient_str,ploidy,i))
        }
      }
    }
  }  
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
saveRDS(df[['ploidy_mix']], sprintf('%s/results/%s/%s_%s_ploidy_mix/%s_%s_ploidy_mix.Rds',user_dir,cohort,cohort,filtering_str,patient_str,filtering_str))
run_lobster(df,patient_str,cohort,filtering_str,user_dir)
      