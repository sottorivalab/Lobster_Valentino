.find_peak_idx <- function(df, bin_width=0.001){
  # find indices and x-values of all peaks found through wavelet-convolution method
  source_python("/Users/madeleine.dale/Documents/TEMULATOR_working/vignettes/find_peaks.py") # change later
  hist_obj <- hist(unlist(df$VAF), plot = FALSE, breaks=seq(from=0, 1, by = bin_width))  
  peak_idx <- find_peaks(hist_obj$counts,bin_width)
  peak_xvals <- hist_obj$breaks[peak_idx]
  return(list(peak_idx = peak_idx, peak_xvals = peak_xvals))
}

.find_minvaf <- function(df, bin_width=0.001){
  # find index and x-value of first peak found through wavelet-convolution method
  peak_idx <- mobster:::.find_peak_idx(df, bin_width=bin_width)$peak_idx
  hist_obj <- hist(unlist(df$VAF), plot = FALSE, breaks=seq(from=0, 1, by = bin_width)) 
  min_vaf <- hist_obj$breaks[peak_idx[1]]
  return(list(min_vaf_idx = peak_idx[[1]], min_vaf = min_vaf))
} 

.log_peak_idx <- function(df){
  # do x log transform, find peaks
  df = df %>% mutate(ln_x = log(VAF))
  peak_idx = mobster:::.find_minvaf(df)$min_vaf_idx
  return(peak_idx)
}

.log_peak_transform <- function(df){
  # do x log transform, find Pareto peak and exclude below, then do min-max scaling
  df = df %>% mutate(ln_x = log(VAF))
  # finding true min, excluding data below
  peaks = mobster:::.find_peak_idx(df)
  peak_idx = peaks$peak_idx
  min_vaf = peaks$peak_xvals[[1]]
  df = df[df$VAF>min_vaf,]
  df = df %>% mutate(ln_x_scaled = mobster:::.custom_MinMaxScale(ln_x,0.0015,0.99)) %>%   
            dplyr::rename(vaf_orig = VAF) %>% mutate(VAF = ln_x_scaled)
  return(list(df = df, peak_idx = peak_idx))
}

peak_select_hist <- function(hist_obj, peak_idx, subset=1:100, bin_width=0.005, untransform_axis = FALSE, df = NULL){
   # return histogram ggplot2 object with bins corresponding to found peaks highlighted
   hist_data <- data.frame(
     midpoints = hist_obj$mids,
     counts = hist_obj$counts,
     probability = hist_obj$density,
     is_peak = ifelse(seq_len(length(hist_obj$counts)) %in% peak_idx, TRUE, FALSE)
   )
   subset_data <- hist_data[subset, ]
   plot_obj <- ggplot(subset_data, aes(x = midpoints, y = probability,fill=is_peak))+
      geom_bar(stat = 'identity', width=bin_width) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      labs(title = sprintf("Peak selection: bin width %f",bin_width), x = "VAF", y = "Prob. density") +
      scale_fill_manual(values = c("steelblue","red"))
   
   if ((untransform_axis)&&(!is.null(df))){
     plot_obj <- plot_obj+
         scale_x_continuous(breaks = seq(0.1,1.0,0.1),  
                            labels = signif(mobster:::.log_untransform(seq(0.1,1.0,0.1),df), digits=2)) }   
  
  return(plot_obj)
}