library(tidyr)
library(dplyr)

config = mobster:::return_config_dict()
user_dir = config$paths$user_dir

cohort = 'VALE'
filtering = 'pass'
directory = sprintf("%s/results/%s/%s_%s_ploidy_mix_new",user_dir,tolower(cohort),tolower(cohort),filtering)
final_directory = paste0(directory,"/final")

files <- list.files(path=final_directory, pattern = "BH.RDS$",full.names = TRUE)
file_strs <- list.files(path=final_directory, pattern = "BH.RDS$",full.names = FALSE)
patient_strs = sub("^([^_]+_[^_]+)_(.*)$", "\\1", file_strs)
ploidy_strs = sub("^([^_]+_[^_]+)_ploidy_(\\d|NA).RDS$", "\\2", file_strs)

all_column_names <- unique(unlist(lapply(files, function(file) names(readRDS(file)))))

# Function to add missing columns to a data frame and fill with NA
add_missing_columns <- function(df,cohort) {
  library(data.table)
  df <- df %>% select(-any_of(c("IS_ESCAPE.x.x", "IS_ESCAPE.y.y")))
  max_t = 3
  missing_columns <- setdiff(1:max_t, as.numeric(names(df[["ploidy"]])))
  setnames(df$"CCF_expected", old = names(df$"CCF_expected"), new = as.character(1:ncol(df$"CCF_expected")))
  setnames(df$"prob_H0", old = names(df$"prob_H0"), new = as.character(1:ncol(df$"prob_H0")))
  setnames(df$"sig_missing", old = names(df$"sig_missing"), new = as.character(1:ncol(df$"sig_missing")))
  for (col in missing_columns) {
    col_name = as.character(col)
    df$"VAF"[[col_name]] <- rep(NA,nrow(df))
    df$"times"[[col_name]] <- rep(NA,nrow(df))
    df$"ploidy"[[col_name]] <- rep(NA,nrow(df))
    df$"ploidy_orig"[[col_name]] <- rep(NA,nrow(df))
    df$"VAF_ln"[[col_name]] <- rep(NA,nrow(df))
    df$"VAF_ln_scaled"[[col_name]] <- rep(NA,nrow(df))
    df$"VAF_orig"[[col_name]] <- rep(NA,nrow(df))
    df$"FILTER"[[col_name]] <- rep(NA,nrow(df))
    df$"DP"[[col_name]] <- rep(NA,nrow(df))
    df$"NV"[[col_name]] <- rep(NA,nrow(df))
    df$"CCF_expected"[[col_name]] <- rep(NA,nrow(df))
    df$"prob_H0"[[col_name]] <- rep(NA,nrow(df))
    df$"sig_missing"[[col_name]] <- rep(NA,nrow(df))
  }
  return(df)
}

# Read each data frame, add missing columns, and store in the list
dataframes <- lapply(seq_along(files), function(i) {
  df <- readRDS(files[i])
  df$IS_ESCAPE.x = NULL
  df$IS_ESCAPE.y = NULL
  df$DP_percent = NULL
  df$percentile = NULL
  df$CCF = NULL
  df$patient = rep(patient_strs[i], nrow(df))
  df$ploidy_fit = rep(ploidy_strs[i], nrow(df))
  df <- add_missing_columns(df,cohort) %>% 
    unnest(cols = VAF, names_sep = '.') %>%
    unnest(cols = times, names_sep = '.') %>%
    unnest(cols = DP, names_sep = '.') %>%
    unnest(cols = NV, names_sep = '.') %>%
    unnest(cols = ploidy, names_sep = '.') %>%
    unnest(cols = ploidy_orig, names_sep = '.') %>%
    unnest(cols = VAF_ln, names_sep = '.') %>%
    unnest(cols = VAF_ln_scaled, names_sep = '.') %>%
    unnest(cols = VAF_orig, names_sep = '.') %>%
    unnest(cols = prob_H0, names_sep = '.') %>%
    unnest(cols = sig_missing, names_sep = '.') %>%
    unnest(cols = CCF_expected, names_sep = '.') %>%
    unnest(cols = FILTER, names_sep = '.')
    
  return(df)
})

merged_df <- do.call(rbind, dataframes)

saveRDS(merged_df, file = paste0(directory,"/merged_final_BH_0.05.Rds"))

# checking that all mutations are accounted for
patient_split = merged_df %>% group_split(patient) 
for (pdf in patient_split){
  print(nrow(pdf))
  if (nrow(pdf)!=unique(pdf$total_muts_all_ploidies)[1]){
    print(unique(pdf$patient)[1])
  }
}

# sample purity table
merged_df$ploidy_fit = sub("^(.*_ploidy_)(.+)_BH.*\\.RDS$", "\\2", merged_df$ploidy_fit)
sample_purity = merged_df %>% 
  dplyr::select(patient, starts_with("VAF_orig"), significance, ploidy_fit,starts_with("DP"), starts_with("times")) %>% 
  dplyr::filter(significance=='TRUNCAL') %>% 
  dplyr::filter(ploidy_fit %in% c('2','NA'))
# filter rows by coverage for better estimate
min_cov = 10
low_cov_rows = sample_purity %>% select(starts_with("DP")) %>% apply(1, function(row) any(row < min_cov, na.rm = TRUE)) 

sample_purity = sample_purity[!low_cov_rows,]
sample_purity <- sample_purity %>% 
  group_by(patient) %>%
  dplyr::reframe(
    across(!starts_with("VAF_orig"), first),  # Keep all columns
    across(starts_with("VAF_orig"), function(x) {2*mean(x, na.rm = TRUE)})
  ) 
colnames(sample_purity) <- gsub("^VAF_orig\\.", "purity.", colnames(sample_purity))
write.csv(sample_purity, file = sprintf("%s/results/%s/%s_sample_purity.csv",user_dir,tolower(cohort),tolower(cohort)))
