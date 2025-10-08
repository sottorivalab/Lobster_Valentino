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
library(yaml)
theme_set(theme_cowplot())

###-------------------------------------------------------
# FUNCTIONS
replace_unfitted_ploidies <- function(ploidy_df,fitted_ploidies){
  
    find_closest_values <- function(orig_ploid, vector) {
      closest_values = c()
      orig_ploid = as.numeric(orig_ploid)
      for (i in orig_ploid){ # for loop because I'm lazy
        lower_vector = vector[i>vector]
        if (length(lower_vector)==0){ # minimum value greater than i
          closest_value = min(vector)
        } else {
          diffs <- i - lower_vector
          closest_value <- lower_vector[which.min(diffs)]
        }
        closest_values = c(closest_values,closest_value)
      }
      return(closest_values)
    }
    
    replace_elements = function(df,need_replacing,replacements){
      for (i in 1:length(need_replacing)){
        df = replace(df,df==need_replacing[i], replacements[i])
      }
      return(df)
    }
  
  # find all original ploidy combinations
  unique_rows <- apply(ploidy_df, 1, paste, collapse = ",") %>% unique()
  unique_df <- data.frame(do.call(rbind, strsplit(unique_rows, ","))) %>% 
    mutate_all(function(x){as.numeric(x)})
  need_replacing = unique_df[sapply(unique_df, function(x){!(x %in% fitted_ploidies)})] %>% unlist()
  replacements = find_closest_values(need_replacing, fitted_ploidies)
  if (length(need_replacing)>0){
    ploidy_df = replace_elements(ploidy_df,need_replacing,replacements) 
  }
  
  return(ploidy_df)
}

mix_ploidy_replace_w_first <- function(df){
  for (i in 2:ncol(df)) {
    df[, i] <- df[, 1]
  }
  return(df)
}

column_varying_ploidies <- function(df){
  
  row_varies <- function(row) {
    length(unique(row,na.rm = FALSE)) > 1
  }
  
  ploidy_df = df$ploidy
  varying_rows <- apply(ploidy_df, 1, row_varies)
  varying_elements = apply(ploidy_df[varying_rows, ], 1, function(row) {
                                paste(row, collapse = ",")})
  df$ploidy_change = rep(NA,nrow(df))
  df$ploidy_change[varying_rows] <- varying_elements
  
  return(df)
}

##################################################
# Main code

config = mobster:::return_config_dict()
user_dir = config$paths$user_dir
from_dataframe = config$from_dataframe

results_dir = sprintf('%s/results',user_dir)
cohort = 'VALE'
imax = 37
filtering = 'pass'

if (from_dataframe==TRUE){
  fit_path = sprintf('%s/%s/%s_fits_from_dataframe/annotated',results_dir,tolower(cohort),tolower(cohort))
} else {
  fit_path = sprintf('%s/%s/%s_fits_pass/annotated',results_dir,tolower(cohort),tolower(cohort))
}
  
for (i in seq(1,imax)){

    prefix <- "transformed_"
    suffix <- ".Rds"
    fit_files <- list.files(fit_path, full.names = TRUE,recursive = FALSE, include.dirs = F)
    fit_files_wo_pref <- list.files(fit_path, full.names = FALSE,recursive = FALSE, include.dirs = F)
    file_strs <- gsub(prefix, "", fit_files_wo_pref)
    file_strs <- gsub(suffix, "", file_strs)
    file_strs <- gsub('wgs_', "", file_strs)
    file_strs <- gsub('_denovo', "", file_strs)
    
    patient_strs = sub("^([^_]+_[^_]+)_(.*)$", "\\1", file_strs)
    ploidy_strs = sub("^([^_]+_[^_]+)_(.*)$", "\\2", file_strs)
    ploidy_strs <- gsub('_pass', "", ploidy_strs)
    ploidy_strs <- gsub('_strandbias', "", ploidy_strs)
    
    patient_str = unique(patient_strs)[i]
    
    ploidy_mix = readRDS(sprintf('%s/%s/%s_pass_ploidy_mix_new/%s_pass_ploidy_mix.Rds',results_dir,tolower(cohort),tolower(cohort),patient_str))
    
    # GATHER FITS FOR DIFF PLOIDIES
    idxs = which(patient_strs==unique(patient_strs)[i])
    file_str_sub = file_strs[idxs]
    ploidy_str_sub = ploidy_strs[idxs] 
    fits = list()
    for (l in 1:length(idxs)){
      fit = readRDS(sprintf('%s/transformed_%s_denovo.Rds',fit_path,file_str_sub[l]))
      ploidy_str = gsub("ploidy_", "", ploidy_str_sub[l])
      fits[[ploidy_str]] = fit
    } 
    fitted_ploidies = names(fits) %>% as.numeric()
      
    if (!is.null(ploidy_mix)){
      ploidy_mix$VAF = ploidy_mix$VAF %>% mutate_all(~ ifelse(is.nan(.), NA, .))
      ploidy_mix <- ploidy_mix[apply(ploidy_mix$VAF, 1, function(x) !all(is.na(unlist(x)))), ]
      print(paste0('nrow(ploidy_mix): ',nrow(ploidy_mix)))
      
      # assume NA in ploidy_mix are diploid - note in new column
      ploidy_mix[['ploidy_orig']]= ploidy_mix[['ploidy']]
      ploidy_mix = column_varying_ploidies(ploidy_mix)   
      ploidy_mix$ploidy = ploidy_mix$ploidy %>% 
          replace(is.na(.), 2)
      
      # if ploidy changes between timepoints, select first ploidy
      ploidy_mix$ploidy = mix_ploidy_replace_w_first(ploidy_mix$ploidy)
      ploidy_mix$ploidy = replace_unfitted_ploidies(ploidy_mix$ploidy, fitted_ploidies) 
      
      # take unique ploidy combinations (should all be same)
      unique_combos <- ploidy_mix$ploidy[,1] %>% unique() %>% unlist()
      # iterate through combinations, annotating using appropriate fit
      for (pli in 1:length(fitted_ploidies)){
        pl = fitted_ploidies[pli] 
      
        if (pl %in% unique_combos){
          df_ploidy = ploidy_mix %>% dplyr::filter(ploidy[[1]]==pl)
          # update to append to unknown if diploid isn't found
          fit = fits[[as.character(pl)]]
          fit$annotated[['ploidy_orig']]= fit$annotated[['ploidy']]
          
          selected_fit = mobster:::select_fit(fit,patient_str,pl,cohort,user_dir)
          new_df = mobster:::Clusters_denovo(selected_fit, df_ploidy)
          
          print(sprintf('ploidy_mix for pl %d: %d',pl,nrow(ploidy_mix)))
          fit$annotated$ploidy_change = rep(NA,nrow(fit$annotated))
          fit$annotated = rbind(fit$annotated,new_df)
        } else {
          fit = fits[[as.character(pl)]]
          fit$annotated[['ploidy_orig']]= fit$annotated[['ploidy']]
          fit$annotated$ploidy_change = rep(NA,nrow(fit$annotated))
        }
        ploidy_str = paste0("ploidy_",pl)
        
        saveRDS(fit, sprintf('%s/%s/%s_pass_ploidy_mix_new/annotated/transformed_%s_%s.RDS',results_dir,tolower(cohort),tolower(cohort),patient_str,ploidy_str))
      }
    } else {
      for (pli in 1:length(fitted_ploidies)){
        pl = fitted_ploidies[pli] 
        # update to append to unknown if diploid isn't found
        fit = fits[[names(fits)[[pli]]]]
        fit$annotated[['ploidy_orig']]= fit$annotated[['ploidy']]
        ploidy_str = paste0("ploidy_",pl)
        fit$annotated$ploidy_change = rep(NA,nrow(fit$annotated))
        
        saveRDS(fit, sprintf('%s/%s/%s_pass_ploidy_mix_new/annotated/transformed_%s_%s.RDS',results_dir,tolower(cohort),tolower(cohort),patient_str,ploidy_str))
      }
    }
#  }
}
