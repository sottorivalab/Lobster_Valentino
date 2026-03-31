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
library(mclust)

#########################################
run_lobster <- function(df,user_dir,fit_str){
  
  xstar = 0.05
  Kvec = seq(3) # max no. of possible beta components 
  
  mfit = df %>% mobster_fit(K=Kvec, samples=1, parallel = TRUE,N_cutoff=1,maxIter = 15, genSA_steps = 3, quantile=TRUE,log_MMscaled = FALSE, untf_convolved = TRUE, xstar = xstar, Gopt = 'DE')
  df$lobster_cluster = mfit$best$data$cluster
  saveRDS(df, file=sprintf('/scratch/madeleine.dale/lobster_paper_fresh/lobster4publishing.VALE/vignettes/fitting_validation/sampling_simulations/dataframes_timeseries/%s.Rds',fit_str))
  saveRDS(mfit, file=sprintf('/scratch/madeleine.dale/lobster_paper_fresh/lobster4publishing.VALE/vignettes/fitting_validation/sampling_simulations/full_dataframes_timeseries/%s.Rds',fit_str))
  
}
###############################################

format_from_TEMULATOR <- function(df){
  
  # for time-series simulation
  df$DP <- df$depth
  df$VAF_orig <- df$VAF
  df$VAF[df$VAF == 0] <- NA
  
  df <- df %>%
    dplyr::mutate(VariantType = "SNP")
  
  return(df)
}

user_dir = ''
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file string)", call.=FALSE)
}

fit_str = args[[1]]

df = readRDS(sprintf("/scratch/madeleine.dale/lobster_paper_fresh/lobster4publishing.VALE/vignettes/fitting_validation/sampling_simulations/dataframes_timeseries/%s.Rds",fit_str))
formatted_df = format_from_TEMULATOR(df)
run_lobster(formatted_df,user_dir,fit_str)
