# simulation sampling from user-defined distribution
library(dplyr)
library(tidyr)

set.seed(123)

# Function to generate Beta values using Mean and Precision
rbeta <- function(n, mean, precision) {
  # Convert Mean/Precision to Shape1/Shape2
  alpha <- mean * precision
  beta  <- (1 - mean) * precision
  
  return(stats::rbeta(n, shape1 = alpha, shape2 = beta))
}

rpareto_sim <- function(n, shape, scale) {
  # n: number of samples
  # shape: alpha 
  # scale: x_min 
  
  u <- runif(n)
  tail_vafs = scale / (1 - u)^(1/shape) # quantile
  tail_vafs = tail_vafs[tail_vafs <= 1.0] # truncate to 1
  
  return(tail_vafs)
}

simulate_clonal_tumour_data <- function(n_clonal=500, n_tail=500, 
                                clonal_prec = 20,
                                pareto_shape = 1, pareto_scale = 0.005) {
  
  # Generate Clonal Cluster (Beta)
  vaf_clonal <- rbeta(n_clonal, 0.5, clonal_prec)
  
  # Generate Tail (Pareto)
  raw_pareto <- rpareto_sim(n_tail, pareto_shape, pareto_scale) 
  vaf_tail   <- raw_pareto[raw_pareto <= 1.0][1:n_tail]
  
  # Combine
  data <- data.frame(
    VAF = c(vaf_clonal, vaf_tail),
    Type = c(rep("Clonal", n_clonal), rep("Tail", n_tail))
  )
  
  return(data)
}

simulate_clonal_subclonal_tumour_data <- function(n_clonal=500, n_subclonal=400, n_tail=500, 
                                        subclonal_ccf = 0.5,
                                        beta_prec = 20,
                                        pareto_shape = 1, pareto_scale = 0.005) {
  
  # Generate Clonal Cluster (Beta)
  vaf_clonal <- rbeta(n_clonal, 0.5, beta_prec) # clonal mean = 0.5
  vaf_subclonal <- rbeta(n_subclonal, 0.5*subclonal_ccf, beta_prec)
   
  # Generate Tail (Pareto)
  raw_pareto <- rpareto_sim(n_tail, pareto_shape, pareto_scale) 
  vaf_tail   <- raw_pareto[raw_pareto <= 1.0][1:n_tail]
  
  # Combine
  data <- data.frame(
    VAF = c(vaf_clonal, vaf_subclonal, vaf_tail),
    Type = c(rep("Clonal", n_clonal), rep("Subclonal", n_subclonal), rep("Tail", n_tail))
  )
  
  return(data)
}

simulate_realistic_reads <- function(ideal_vafs, purity, depth_mean = 300, depth_dispersion=25) {
  
  n <- nrow(ideal_vafs)
  
  # Scale by purity 
  true_vafs <- ideal_vafs$VAF * purity
  
  # Simulate over-dispersed sequencing Depth (negative binomial)
  # mu = Average coverage
  # size = Inverse dispersion (smaller number = wider spread)
  depths <- rnbinom(n, mu = depth_mean, size = depth_dispersion)
  #print(hist(depths, breaks=15))
  
  # Simulate Reads (Binomial)
  alt_reads <- rbinom(n, size = depths, prob = true_vafs)
  
  df = data.frame(
    depth = depths,
    alt = alt_reads,
    VAF = alt_reads / depths,
    true_vafs = true_vafs,
    Type = ideal_vafs$Type
  )
  
  #df = df %>% dplyr::filter(VAF!=0)
  
  return(df)
}

combine_to_nested <- function(df_list) {
  
  combined_df <- data.frame(Type = df_list[[1]]$Type)
  
  cols_to_nest <- c("depth", "alt", "VAF", "true_vafs")
  
  for (col in cols_to_nest) {
    
    temp_nested_df <- data.frame(matrix(ncol = length(df_list), nrow = nrow(combined_df)))
    
    names(temp_nested_df) <- as.character(seq_along(df_list))
    
    for (i in seq_along(df_list)) {
      temp_nested_df[[as.character(i)]] <- df_list[[i]][[col]]
    }
    
    combined_df[[col]] <- temp_nested_df
  }
  
  return(combined_df)
}

output_dir = '/scratch/madeleine.dale/lobster_paper_fresh/lobster4publishing.VALE/vignettes/fitting_validation/sampling_simulations'

true_clonal_vafs = simulate_clonal_tumour_data()
subclonal_ccfs = c(0, 0.2, 0.35, 0.5, 0.65) # 0 = clonal simulation
purity_1 = c(0.5, 0.25, 0.1) # timepoint 1
purity_2 = c(0.5, 0.25, 0.1, 0.05)
purity_3 = c(0.5, 0.25, 0.1)

all_combos <- crossing(subclonal_ccfs = subclonal_ccfs, purity_1 = purity_1, purity_2 = purity_2, purity_3 = purity_3)

results_list <- pmap(all_combos, function(subclonal_ccfs, purity_1, purity_2, purity_3) {
  
  print(sprintf("Sampling params: subclonal_ccf: %f, purity_1: %f, purity_2: %f, purity_3: %f", subclonal_ccfs, purity_1, purity_2, purity_3))
  
  if (subclonal_ccfs==0){
    true_vafs = simulate_clonal_tumour_data()
  } else {
    true_vafs = simulate_clonal_subclonal_tumour_data(subclonal_ccf = subclonal_ccfs)
  }
   
  timepoint_1 = simulate_realistic_reads(true_vafs, purity = purity_1)
  timepoint_2 = simulate_realistic_reads(true_vafs, purity = purity_2)
  timepoint_3 = simulate_realistic_reads(true_vafs, purity = purity_3)
  
  df = combine_to_nested(list(timepoint_1, timepoint_2, timepoint_3))
  df = df[rowSums(df$VAF) > 0, ]
  
  saveRDS(df, file = sprintf("/scratch/madeleine.dale/lobster_paper_fresh/lobster4publishing.VALE/vignettes/fitting_validation/sampling_simulations/dataframes_timeseries/subclonal_ccf_%f_purity1_%f_purity2_%f_purity3_%f.Rds", subclonal_ccfs, purity_1, purity_2, purity_3))
  return(df)
})

for (subclone_ccf in c(0.2, 0.35, 0.5, 0.65)){
  
  for (purity in c(1.0, 0.5, 0.25, 0.1)){
    true_subclonal_vafs = simulate_clonal_subclonal_tumour_data(subclonal_ccf = subclone_ccf)
    hist(true_subclonal_vafs$VAF, breaks = 50, col = "skyblue", main = "Ideal subclonal VAF distribution")
    
    sequenced_subclonal_vafs = simulate_realistic_reads(true_subclonal_vafs, purity = purity)
    #pdf(sprintf("%s/histograms/sequenced_subclone_ccf_%f_purity_%f.pdf",output_dir,subclone_ccf,purity), width = 8, height = 6)
    #hist(sequenced_subclonal_vafs$VAF, breaks = 50, col = "skyblue", main = "Sequenced subclonal VAF distribution")
    #dev.off()
    
    saveRDS(sequenced_subclonal_vafs, file=sprintf('%s/dataframes_timeseries/subclone_ccf_%f_purity_%f.Rds',output_dir,subclone_ccf,purity))
  }
}

for (purity in c(1.0, 0.5, 0.25, 0.1)){
  
    hist(true_clonal_vafs$VAF, breaks = 50, col = "skyblue", main = "Ideal clonal VAF distribution")
    
    sequenced_clonal_vafs = simulate_realistic_reads(true_clonal_vafs, purity = purity)
    #pdf(sprintf("%s/histograms/sequenced_clone_purity_%f.pdf",output_dir,purity), width = 8, height = 6)
    #hist(sequenced_clonal_vafs$VAF, breaks = 50, col = "skyblue", main = "Sequenced clonal VAF distribution")
    #dev.off()
    sequenced_clonal_vafs = simulate_realistic_reads_timeseries
    
    saveRDS(sequenced_clonal_vafs, file=sprintf('%s/dataframes_timeseries/clone_purity_%f.Rds',output_dir,purity))
    
}
