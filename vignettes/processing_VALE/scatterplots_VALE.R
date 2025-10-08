library(cowplot)
library(ggplot2)
library(plyr)
library(dplyr)
library(doParallel)
library(mobster)
library(R.utils)
library(rtern)
library(fields)
library(RColorBrewer)
library(ggrepel)
library(tidyr)
library(purrr)
theme_set(theme_cowplot())
library(ggpubr)
library(gridExtra)
library(yaml)

###--------------------------------------------------------------------------------------------------------------------------------------
# FUNCTIONS
plot_paper = function(x, timepoint_labels,
                      labeled_colours,
                      cutoff_assignment = 0,
                      annotation_extras = NULL,
                      secondary_axis = NULL,
                      MM_scale_list=NULL, 
                      log_MMscaled=FALSE,
                      transform_axis_back=FALSE,
                      highlight_NA = FALSE,
                      annotate_drivers = FALSE,
                      annotate_escape = FALSE,
                      xlims = NULL,
                      ylims = NULL,
                      ...)
{
  MM_scale_a = MM_scale_list$MM_scale_a
  MM_scale_b = MM_scale_list$MM_scale_b
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Plotting variables
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  binwidth = 0.01
  domain = seq(0, 1, binwidth)
  
  hist_pl = lapply(as.list(1:ncol(x$times)),FUN = function(t){individual_plot_paper(t = t, timepoint_label = timepoint_labels[t], x = x, 
                                                                      labeled_colours = labeled_colours,
                                                                      domain = domain,
                                                                      binwidth = binwidth,
                                                                      cutoff_assignment = cutoff_assignment,
                                                                      beta_colors = beta_colors,
                                                                      tail_color = tail_color,
                                                                      annotation_extras = annotation_extras,
                                                                      secondary_axis = secondary_axis,
                                                                      MM_scale_list = MM_scale_list, 
                                                                      log_MMscaled = log_MMscaled,
                                                                      transform_axis_back = transform_axis_back,
                                                                      highlight_NA = highlight_NA,
                                                                      annotate_drivers = annotate_drivers,
                                                                      annotate_escape = annotate_escape,
                                                                      xlims = xlims,
                                                                      ylims = ylims
  )
  })
  
  return(hist_pl)
}

individual_plot_paper <- function(t, timepoint_label, x, labeled_colours, domain, binwidth,
                            cutoff_assignment = 0,
                            beta_colors = RColorBrewer::brewer.pal(n = 9, 'Set1'),
                            na_colors = RColorBrewer::brewer.pal(n = 8, 'Set2'),
                            driver_colors = RColorBrewer::brewer.pal(n = 8, 'Set3'),
                            tail_color = "gainsboro",
                            annotation_extras = NULL,
                            secondary_axis = NULL,
                            MM_scale_list=NULL, 
                            log_MMscaled=FALSE,
                            transform_axis_back=FALSE,
                            highlight_NA = FALSE,
                            annotate_drivers = FALSE,
                            annotate_escape = FALSE,
                            xlims = NULL,
                            ylims = NULL
)
{
  MM_scale_a = MM_scale_list$MM_scale_a
  MM_scale_b = MM_scale_list$MM_scale_b
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Plotting data
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Clustering assignments and clusters
  
  x$VAF = if(log_MMscaled==TRUE) x$VAF_ln_scaled else x$VAF_orig 
  plot_data = x %>% dplyr::arrange(significance)
  
  clusters = sort(unique(plot_data$significance), na.last = TRUE)
 
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Plotting aesthetics
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  # Plot title
  plot_title = sprintf('VAF %s', timepoint_label)
  plot_subtitle = paste0("n = ", nrow(plot_data[!is.na(plot_data$VAF[[as.character(t)]]), ]))
  
  
  total_count <- nrow(plot_data[complete.cases(plot_data$VAF[[as.character(t)]]), ])
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Main ggplot object is the histogram
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  hist_pl = ggplot(as_tibble(plot_data),
                   aes(x = VAF[[as.character(t)]],
                       fill = factor(significance),
                       y = ..count.. / total_count 
                   )) +
    scale_fill_manual(values = labeled_colours) + 
    geom_histogram(
      alpha = .6,
      color = NA,
      position = 'stack',#'identity',
      binwidth = binwidth
    ) + 
    guides(color = FALSE) +
    labs(
      title = bquote(.(plot_title)),
      x = "Observed VAF",
      y = "Density"
    ) +
    annotate("text", x = Inf, y = Inf, label = plot_subtitle, 
               hjust = 1, vjust = 1, size = 5) + 
    theme_minimal() + 
    coord_cartesian(clip = 'off') + 
    labs(fill = "Clusters") 
  
    my_data = as_tibble(plot_data) %>% filter(IS_DRIVER==TRUE)
    
    if ((all(c('GENE','IS_DRIVER') %in% names(as_tibble(plot_data))))&
        annotate_drivers){
      hist_pl = hist_pl + 
        geom_text_repel(
          data = as_tibble(plot_data) %>% filter(IS_DRIVER==TRUE),
          aes(label = GENE, x = VAF[[as.character(t)]], 
              y = rep(0,nrow(as_tibble(plot_data) %>% filter(IS_DRIVER==TRUE))),
              fontface = "bold"), 
          inherit.aes = FALSE,
          max.overlaps = 100,
          size = 2.7, # text size
          force=5,
          nudge_y = 0.01,
          force_pull = 1,
          segment.color = "grey40", # line colour
          segment.size = 0.1,
          angle = 90
        )  
    }
  
    if ((all(c('GENE','IS_ESCAPE') %in% names(as_tibble(plot_data))))&
        annotate_escape){
      hist_pl = hist_pl + 
        geom_text_repel(
          data = as_tibble(plot_data) %>% filter(IS_ESCAPE==TRUE),
          aes(label = GENE, x = VAF[[as.character(t)]], 
              y = rep(0,nrow(as_tibble(plot_data) %>% filter(IS_ESCAPE==TRUE))),
              fontface = "bold"), 
          inherit.aes = FALSE,
          max.overlaps = 100,
          size = 3, # text size
          force=5,
          nudge_y = 0.01,
          force_pull = 1,
          segment.color = "grey40", # line colour
          segment.size = 0.1,
          angle = 90
        )  
    }
  
   legend = get_legend(hist_pl)
   hist_pl = hist_pl + theme(legend.position="none")
  
  if (transform_axis_back){
    my_seq = seq(0.1,1,0.15)
      hist_pl = hist_pl +
        scale_x_continuous(breaks = my_seq,  
                           labels = signif(mobster:::.log_untransform(my_seq,x), digits=2)) 
  }  
   
   
  
  return(list(hist = hist_pl, legend = legend))
}

tp_string <- function(vector) {
  replacement_map <- c("BL", "BR", "PD")
  for (i in 1:length(replacement_map)){
    vector[vector == i] <- replacement_map[i]
  }
  return(vector)
}

rename_NA_clusters_vale <- function(df, offset){
  
  find_values_not_in_column <- function(val, vector) {
    appearing_values = c()
    if (!is.na(val)){
      for (i in 1:length(vector)){
        if (grepl(vector[i], val)){ appearing_values = c(appearing_values,vector[i]) }
      }
      non_appearing_values = setdiff(vector, appearing_values)
      # overwrite if clonal
      if (grepl('C[1-9]|Tail', val)){ non_appearing_values = val }
    } else { non_appearing_values <- NA }
    
    return(non_appearing_values)
  }
  
  select_clinical_clusters <- function(df){
    NA_elements = df %>% filter(is.na(NA_complement))
    only_appear <- df[!grepl("C[1-9]|Tail", df$NA_complement), ] %>% filter(!is.na(NA_complement),lengths(NA_complement) == 1) %>% mutate(NA_complement = paste0("Only in ", NA_complement)) 
    clonal = df[grepl("C[1-9]|Tail", df$NA_complement), ] 
    
    # what is the last timepoint existing?
    times = unique(df$times)
    possible_times = c('BL','BR','PD') 
    time_strings = possible_times[unlist(times)]
    last_timepoint = time_strings[length(time_strings)]
    existing_combos = list(time_strings)
    well_defined = rbind(only_appear, clonal, NA_elements)
    
    # mutation that become NA in intermediate timepoints but appear at last timepoint need to be separated and will be clonal  
    combos = list(
                  unique(list('BL',last_timepoint)), 
                  unique(list('BR',last_timepoint)),
                  unique(list('BL','BR'))
                 )
    for (combo in combos){
      if (length(combo)>1) { # protect against 'Only in'
        combo_df_idx = df$NA_complement %>% sapply(function(x) identical(x, unlist(combo))) 
        is_in_list <- any(sapply(existing_combos, identical, combo))
        if (!is_in_list){
          combo_df = df[combo_df_idx,] %>% mutate(NA_complement = sprintf('In %s (missing elsewhere)', paste(combo, collapse = ", ")))
          existing_combos = append(existing_combos, list(combo))
          well_defined = rbind(well_defined, combo_df)
        } 
      }
    } 
    
    not_appearing = anti_join(df, well_defined, by=c('POS','CHROM','REF','ALT')) %>% mutate(NA_complement = 'Not relevant')
    reconstructed = rbind(well_defined, not_appearing) %>% dplyr::rename('Clinical cluster' = NA_complement)
    
    return(reconstructed) 
  }
  
  times = unique(df$times)
  possible_times = c('BL','BR','PD') 
  time_strings = possible_times[unlist(times)]
  
  # first convert to actual numerical timepoint
  for (i in seq(1,length(times))){
    df$clusterNA = gsub(sprintf('_%d',i+offset), sprintf('_%d',times[[as.character(i)]]), df$clusterNA)
  }
  
  # then to timepoint name
  for (i in seq(1,max(times))){
    df$clusterNA = gsub(sprintf('_%d',i), sprintf('_%s',possible_times[i]), df$clusterNA)
  }
  
  # find complement of NAs
  df$NA_complement <- sapply(df$clusterNA, FUN=function(x){find_values_not_in_column(x,vector = time_strings)})
  
  df = select_clinical_clusters(df)
  
  return(df)
}

# processing clinical significance table
missoni_clinical_significance = function(df, ploidy, sample, cohort, user_dir){
  
  times = unique(df$times)
  possible_times = c('BL','BR','PD') 
  time_strings = possible_times[unlist(times)]
  last_timepoint = time_strings[length(time_strings)]
  
  sig = read.csv(sprintf('%s/%s/%s_clinical_significance.csv',user_dir,tolower(cohort),tolower(cohort)), header = TRUE, sep=';')
  sig = sig %>% dplyr::rename(SELECTED = POSITIVELY.SELECTED..RESISTANCE.)
  sig$SELECTED <- sapply(sig$SELECTED, function(x) { if (!is.na(x)){trimws(unlist(strsplit(x, split = '\\+')))} else {NA}} ) 
  sig$CLONAL <- sapply(sig$CLONAL, function(x) { if (!is.na(x)){trimws(unlist(strsplit(x, split = '\\+')))} else {NA}})
  sig$SUBCLONAL_TAIL <- sapply(sig$SUBCLONAL, function(x) { if (!is.na(x)){trimws(unlist(strsplit(x, split = '\\+')))} else {NA}})
  sig$SAMPLE = gsub(' ', '_', sig$SAMPLE)
  
  ploidy <- ifelse(is.na(ploidy), "NA", ploidy) %>% as.character()
  sig$PLOIDY = tidyr::replace_na(as.character(sig$PLOIDY),"NA") 
  df$'Clinical cluster' = df$'Clinical cluster' %>% trimws()
  
  df$sig_cluster = df$`Clinical cluster` %>% trimws()
  df$significance = rep(NA,nrow(df))
  
  # select groups that are in TRUNCAL category
  clonal_groups = sig %>% 
    filter(PLOIDY == ploidy, SAMPLE == sample) %>%
    dplyr::select(CLONAL) %>% unlist()
  
  subclonal_groups = sig %>% 
    filter(PLOIDY == ploidy, SAMPLE == sample) %>%
    dplyr::select(SUBCLONAL_TAIL) %>% unlist()
  
  selected_groups = sig %>% 
    filter(PLOIDY == ploidy, SAMPLE == sample) %>%
    dplyr::select(SELECTED) %>% unlist()
  
  # filling in 'automatic' categories
  # missoni 65 (only one with one timepoint) is at PD
  df$significance[(df$'Clinical cluster'=='In BL, PD (missing elsewhere)')&(df$cluster_orig %in% clonal_groups)] = 'TRUNCAL'
  df$significance[(df$'Clinical cluster'=='In BL, PD (missing elsewhere)')&(df$cluster_orig %in% subclonal_groups)] = 'SUBCLONAL/TAIL'
  df$significance[(df$'Clinical cluster'=='In BL, PD (missing elsewhere)')&(df$cluster_orig %in% selected_groups)] = 'EXPANDED'
  
  rows_to_fill = which(df$'Clinical cluster'=='In BL, BR (missing elsewhere)')
  if ('PD' %in% time_strings){
    df$significance[rows_to_fill] = 'CONTRACTED' 
  } else {
    df$significance[(df$'Clinical cluster'=='In BL, BR (missing elsewhere)')&(df$cluster_orig %in% clonal_groups)] = 'TRUNCAL'
    df$significance[(df$'Clinical cluster'=='In BL, BR (missing elsewhere)')&(df$cluster_orig %in% subclonal_groups)] = 'SUBCLONAL/TAIL'
    df$significance[(df$'Clinical cluster'=='In BL, BR (missing elsewhere)')&(df$cluster_orig %in% selected_groups)] = 'EXPANDED'
  }
  
  rows_to_fill = which(df$'Clinical cluster'=='In BR, PD (missing elsewhere)')
  if ('BL' %in% time_strings){
    df$significance[rows_to_fill] = 'EXPANDED' 
  } else {
    df$significance[(df$'Clinical cluster'=='In BR, PD (missing elsewhere)')&(df$cluster_orig %in% clonal_groups)] = 'TRUNCAL'
    df$significance[(df$'Clinical cluster'=='In BR, PD (missing elsewhere)')&(df$cluster_orig %in% subclonal_groups)] = 'SUBCLONAL/TAIL'
    df$significance[(df$'Clinical cluster'=='In BR, PD (missing elsewhere)')&(df$cluster_orig %in% selected_groups)] = 'EXPANDED'
  }
  
  rows_to_fill = which(df$'Clinical cluster'=='Only in BL')
  df$significance[rows_to_fill] = 'CONTRACTED'  # BL timepoint never appears on its own in MISSONI (check VALE!)
  
   # BR timepoint never appears on its own in MISSONI (check VALE!)
  rows_to_fill = which(df$'Clinical cluster'=='Only in BR')
   if ('BL' %in% time_strings) {
    # PD missing  
    df$significance[rows_to_fill] = 'EXPANDED' 
  } else if ('PD' %in% time_strings) { 
    # BL missing
    df$significance[rows_to_fill] = 'CONTRACTED'  # could also be labelled intermediate
  }
  # overwrite if both inside (most common)
  if (all(c('BL','PD') %in% time_strings)){
    df$significance[rows_to_fill] = 'TRANSIENT' 
  }
  
  rows_to_fill = which(df$'Clinical cluster'=='Only in PD')
  if (any('BL' %in% time_strings,'BR' %in% time_strings)){
    df$significance[rows_to_fill] = 'EXPANDED' 
  } else { # unless only PD timepoint -> clonal
    df$significance[(df$'Clinical cluster'=='Only in PD')&(df$cluster_orig %in% clonal_groups)] = 'TRUNCAL'
    df$significance[(df$'Clinical cluster'=='Only in PD')&(df$cluster_orig %in% subclonal_groups)] = 'SUBCLONAL/TAIL'
    df$significance[(df$'Clinical cluster'=='Only in PD')&(df$cluster_orig %in% selected_groups)] = 'EXPANDED'
  }
  
  # no missing timepoint - mutation belongs to clone/subclone 
  df$significance[df$sig_cluster %in% clonal_groups] = 'TRUNCAL' 
  df$significance[df$sig_cluster %in% subclonal_groups] = 'SUBCLONAL/TAIL' 
  df$significance[df$sig_cluster %in% selected_groups] = 'EXPANDED' 
  
  if (!(sample %in% sig$SAMPLE)){df$significance = rep(NA,nrow(df))}
  
  df$sig_cluster = NULL
  
  return(df)
}

find_values_not_in_column <- function(val, vector) {
  appearing_values = c()
  if (!is.na(val)){
    for (i in 1:length(vector)){
      if (grepl(vector[i], val)){ appearing_values = c(appearing_values,vector[i]) }
    }
    non_appearing_values = setdiff(vector, appearing_values)
    # overwrite if clonal
    if (grepl('C[1-9]|Tail', val)){ non_appearing_values = val }
  } else { non_appearing_values <- NA }
  
  return(non_appearing_values)
}
  
# return plot_grid object to wrap scatterplots if too long
wrapped_plots <- function(pdf_plots,legend){
  if (length(pdf_plots)>3){
    my_nrow = 2
    my_ncol = round(length(pdf_plots)/2)
  } else {
    my_nrow = 1
    my_ncol=length(pdf_plots)
  }
  
  wrapped = plot_grid(plot_grid(plotlist=pdf_plots,nrow=my_nrow,align = 'hv',axis = 'tblr', rel_widths = rep(1,length(pdf_plots)), rel_heights = rep(1,length(pdf_plots))), legend,ncol=2,align='hv',rel_widths = c(3.5, 1.0))
  
  return(list(wrapped = wrapped, ncol = my_ncol, nrow = my_nrow))
}
    
to_sentence_case <- function(text) {
  text <- tolower(text)
  text <- paste0(toupper(substr(text, 1, 1)), substr(text, 2, nchar(text)))
  return(text)
}

extract_elements <- function(lst, element_name) {
  extracted_elements <- map(lst, element_name) %>% compact()
  return(extracted_elements)
}

hypothesis_testing <- function(idxs, ploidy_str_sub, fit_path, file_str_sub){
  
    fits = list()
    MM_scales = list()
    
    # hypothesis testing
    for (l in 1:length(idxs)){
      ploidy_str = ploidy_str_sub[l]
      pl = gsub("ploidy_", "", ploidy_str)
      pl = gsub("_annotated_mix", "", pl) %>% as.numeric()
      
      # select appropriate fit
      fit = readRDS(sprintf('%s/transformed_%s.RDS',fit_path,file_str_sub[l])) 
      fit$annotated$VAF <- fit$annotated$VAF %>% mutate_all(~ ifelse(is.nan(.), NA, .)) # replace NaN with NA
      
      # calculate log transformation (separate from hypothesis testing - which uses only VAF_orig)
      MM_scale  = mobster:::get_a_b(log(fit$annotated$VAF))
      fit$MM_scale$MM_scale_a = MM_scale$a 
      fit$MM_scale$MM_scale_b = MM_scale$b 
      MM_scales[[l]] = mobster:::get_a_b(log(fit$best$data$VAF))
      fit$annotated = mobster:::.log_transform_a_b(fit$annotated, 
                                                   a = fit$MM_scale$MM_scale_a, 
                                                   b = fit$MM_scale$MM_scale_b)
      fit$annotated = fit$annotated %>% 
        dplyr::rename(cluster = cluster_denovo) %>% 
        mutate(ploidy_fit = rep(pl,nrow(fit$annotated)))
      fit$annotated = fit$annotated %>% mobster:::MC_GEL_driver_labelling_posthoc()
      
      # calculate purities
      if (ploidy_str %in% c('ploidy_2', 'ploidy_NA')){
        # finding purities
        purities = mobster:::calc_VAF_clone(fit) 
        purities = purities*2 # x2 for sample purity from diploid
      } 
      fits[[l]] = fit
    }
    
    for (l in 1:length(idxs)){
      fit = fits[[l]]
      
      # treat ploidy_NA as ploidy_2
      ploidy_str = ploidy_str_sub[l]
      if (ploidy_str=='ploidy_NA'){ ploidy_str = 'ploidy_2' }
      pl = gsub("ploidy_", "", ploidy_str)
      pl = gsub("_annotated_mix", "", pl) %>% as.numeric()
      
      times = unique(fit$annotated$times)
      # for checking total number of mutations later
      fit$annotated$total_muts_all_ploidies = rep(fit$mutations_over_all_ploidies,nrow(fit$annotated))
      
      # calculate expected CCF according to null hypothesis
      fit$annotated = mobster:::calc_percentile(fit$annotated)
      fit$annotated$CCF_expected = mobster:::calc_CCF_H0(fit)
      
      # replace NAs with zeros
      fit$annotated$DP[is.na(fit$annotated$DP)] = 0 
      fit$annotated$NV[is.na(fit$annotated$NV)] = 0 
      
      # manipulate to correct dimensions
      purity_df = data.frame(matrix(rep(purities, nrow(fit$annotated)), nrow = nrow(fit$annotated), byrow = TRUE))
      pl_df = data.frame(matrix(rep(pl, nrow(fit$annotated)*ncol(fit$annotated$times)), nrow = nrow(fit$annotated), byrow = FALSE))
      
      # calculate p values
      func_list = mobster:::binomial_mixture_probdist(pl) # return appropriate mixture of binomial functions
      # applies function for evaluating prob of zero count under clonality elementwise to dataframes
      fit$annotated$prob_H0 = as.data.frame(Map(function(depth,purity,CCF,pl){  
                      mapply(function(col1,col2,col3,col4){
                        mobster:::safe_eval_funclist_elementwise(func_list,col1,col2,col3,col4)},
                        depth,purity,CCF,pl,SIMPLIFY=TRUE)}, 
                      fit$annotated$DP, purity_df, fit$annotated$CCF_expected, pl_df))
      
      if (l==1){
        fits_all_ploidies = fit$annotated
      } else {
        fits_all_ploidies = rbind(fits_all_ploidies,fit$annotated)
      }
      
      fits[[l]] = fit$annotated
    }
    
    stopifnot((nrow(fits_all_ploidies) == fit$annotated$total_muts_all_ploidies)) # check no lost mutations
    
    # do BH procedure over all ploidies
    fits_all_ploidies$sig_missing = mobster:::sig_BH_proc(fits_all_ploidies$VAF_orig,fits_all_ploidies$prob_H0)
    
    return(list(fits_all_ploidies=fits_all_ploidies, MM_scales=MM_scales, times=times))
}

###-------------------------------------------------------------------------------------------------------------------------------
# MAIN CODE

# global vars ------------------------------------------------

config = mobster:::return_config_dict()
user_dir = config$paths$user_dir

cohort <- 'VALE'
transformed <- FALSE
coverage_filtering <- FALSE
min_cov <- 10

# Determine filtering mode and string
filtering <- 'pass' 
filtering_str <- 'pass' 

# Define paths based on filtering mode
base_dir <- sprintf('%s/results/%s', user_dir, tolower(cohort))
path <- sprintf('%s/%s_scatterplots_%s', base_dir, tolower(cohort), filtering)
fit_path <- sprintf('%s/%s_%s_ploidy_mix_new/annotated', base_dir, tolower(cohort), filtering)
final_fit_path <- sprintf('%s/%s_%s_ploidy_mix_new/final', base_dir, tolower(cohort), filtering)

output_path <- path
prefix <- "transformed_"
suffix <- ".RDS"

# Get and clean file names
fit_files <- list.files(fit_path, full.names = TRUE, recursive = FALSE, include.dirs = FALSE)
fit_files_wo_pref <- list.files(fit_path, full.names = FALSE, recursive = FALSE, include.dirs = FALSE)
file_strs <- gsub(prefix, "", fit_files_wo_pref)
file_strs <- gsub(suffix, "", file_strs)

patient_strs <- sub("^([^_]+_[^_]+)_(.*)$", "\\1", file_strs)
ploidy_strs <- sub("^([^_]+_[^_]+)_(.*)$", "\\2", file_strs)
#-----------------------------------------------------------

# main loop
for (sc in c('scaled', 'unscaled')) { # both transformed and untransformed VAF plots
  for (patient_str in unique(patient_strs)) {
    
    print(sprintf('patient: %s', patient_str))
    
    # select patient's files and ploidies
    idxs <- which(patient_strs == patient_str)
    file_str_sub <- file_strs[idxs]
    ploidy_str_sub <- ploidy_strs[idxs]
    
    # reorder so that, if diploid fit exists, it is at start
    if ('ploidy_2' %in% ploidy_str_sub) {
      ploidy_strs_ordered <- c('ploidy_2', ploidy_str_sub[!ploidy_str_sub %in% 'ploidy_2'])
      permutation <- match(ploidy_strs_ordered, ploidy_str_sub)
      file_str_sub <- file_str_sub[permutation]
      ploidy_str_sub <- ploidy_strs_ordered
    }
    
    # calculate variables via hypothesis testing function
    return_list <- hypothesis_testing(idxs, ploidy_str_sub, fit_path, file_str_sub)
    fits_all_ploidies <- return_list$fits_all_ploidies
    MM_scales <- return_list$MM_scales
    times <- return_list$times
    
    # Define colors once
    avail_colours <- c(RColorBrewer::brewer.pal(5, "Set1")[3],
                       RColorBrewer::brewer.pal(5, "Set3")[5],
                       RColorBrewer::brewer.pal(5, "Set3")[4],
                       RColorBrewer::brewer.pal(11, "Pastel1")[4],
                       RColorBrewer::brewer.pal(11, "Pastel1")[5],
                       RColorBrewer::brewer.pal(8, "Accent"),
                       RColorBrewer::brewer.pal(8, "Dark2"),
                       RColorBrewer::brewer.pal(12, "Set3"),
                       RColorBrewer::brewer.pal(8, "Set2"))
    
    my_levels3 <- c('TRUNCAL', 'SUBCLONAL/TAIL', 'EXPANDED', 'TRANSIENT', 'CONTRACTED')
    my_colours3 <- setNames(avail_colours[1:5], my_levels3)
    
    pdf_plots_all9 <- list()
    pdf_plots_all10 <- list()
    
    for (l in 1:length(idxs)) {
      ploidy_str <- ploidy_str_sub[l]
      pl <- gsub("ploidy_", "", ploidy_str) %>%
        gsub("_annotated_mix", "", .) %>%
        as.numeric()
      
      # Select appropriate fit
      if (!is.na(pl)) {
        fit <- fits_all_ploidies %>% filter(ploidy_fit == pl)
      } else {
        fit <- fits_all_ploidies %>% filter(is.na(ploidy_fit))
      }
      
      # Apply cohort-specific posthoc functions
      possible_times <- c('BL', 'BR', 'PD')
      time_strings <- possible_times[unlist(times)]
      offset <- which(possible_times == time_strings[1]) - 1
      fit <- mobster:::posthoc_NA_clusters(fit, offset, mask = TRUE)
      fit <- rename_NA_clusters_vale(fit, offset)
      fit <- missoni_clinical_significance(fit, pl, patient_str, cohort, user_dir)
      if (patient_str == 'VALE_003-006') {
        fit[fit$GENE == 'KRAS', ]$significance <- 'EXPANDED'
      }
      
      # Save RDS file
      fit %>%
        dplyr::select(-IS_DRIVER, -IS_ESCAPE) %>%
        saveRDS(file = sprintf('%s/%s_BH.RDS', final_fit_path, file_str_sub[l]))
      
      # Exclude rows with all NA VAFs
      na_all_idx <- which(rowSums(is.na(fit$VAF_orig)) == ncol(fit$VAF_orig))
      if (length(na_all_idx) > 0) {
        fit <- fit[-na_all_idx, ]
      }
      
      # Coverage filtering
      if (coverage_filtering == TRUE) {
        fit$VAF[fit$DP < min_cov] <- NA
        fit$VAF_orig[fit$DP < min_cov] <- NA
        fit$VAF_ln_scaled[fit$DP < min_cov] <- NA
      }
      
      # Generate plots for drivers and escape genes
      for (plot_type in c('drivers', 'escape')) {
        
        # Determine annotation based on plot type
        annotate_drivers <- plot_type == 'drivers'
        annotate_escape <- plot_type == 'escape'
        
        # Get plot objects
        hist_plot_objects <- list(
          plot_paper(fit, time_strings, my_colours3, log_MMscaled = (sc == 'scaled'), transform_axis_back = (sc == 'scaled'), MM_scale_list = MM_scales[[l]], highlight_NA = FALSE),
          plot_paper(fit, time_strings, my_colours3, log_MMscaled = (sc == 'scaled'), transform_axis_back = (sc == 'scaled'), MM_scale_list = MM_scales[[l]], highlight_NA = FALSE, annotate_drivers = annotate_drivers),
          plot_paper(fit, time_strings, my_colours3, log_MMscaled = (sc == 'scaled'), transform_axis_back = (sc == 'scaled'), MM_scale_list = MM_scales[[l]], highlight_NA = FALSE, annotate_escape = annotate_escape)
        )
        
        tables <- mobster:::get_driver_escape_table(fit, cohort)
        
        # Define plot titles
        title <- ggdraw() +
          draw_label(
            paste0(to_sentence_case(gsub("_", " fit: ", ploidy_str)), if (coverage_filtering) " (coverage filtered!)" else ""),
            fontface = 'bold', x = 0, hjust = 0, size = 20
          ) + theme(plot.margin = margin(0, 0, 0, 7))
        
        subtitle <- ggdraw() +
          draw_label(
            sprintf('%s genes labelled', ifelse(plot_type == 'drivers', 'Driver', 'Escape')),
            x = 0, hjust = 0, size = 15
          ) + theme(plot.margin = margin(0, 0, 7, 7))
        
        hist_plot_object <- hist_plot_objects[[ifelse(annotate_drivers, 2, 3)]]
        hist_plots <- extract_elements(hist_plot_object, "hist")
        hist_plot_object <- plot_grid(plotlist = hist_plots, ncol = length(hist_plots), nrow = 1)
        
        # Generate and save plots based on number of timepoints
        avail_timepoints <- unlist(fit$times[1, ]) %>% as.vector()
        
        if (length(avail_timepoints) < 2) {
          intermed <- plot_grid(
            hist_plot_object,
            plot_grid(extract_elements(hist_plot_objects[[1]], "legend")[[1]], tables[[plot_type]], nrow = 2, rel_heights = c(0.3, 0.7), align = "v"),
            ncol = 2, rel_widths = c(3.5, 0.7)
          )
          arranged <- plot_grid(title, intermed, ncol = 1, nrow = 2, rel_heights = c(0.05, 1))
        } else {
          # Scatter plot logic
          my_combs <- combn(avail_timepoints, 2)
          time_combos <- if (!is.null(ncol(my_combs))) split(my_combs, rep(1:ncol(my_combs), each = nrow(my_combs))) %>% unname() else list(my_combs)
          axis_titles <- lapply(time_combos, tp_string) #, cohort = cohort
          
          # Prepare data frame for plotting
          my_data <- if (sc == "scaled") fit$VAF_ln_scaled else fit$VAF_orig
          BL_MR_PD <- as.data.frame(cbind(fit$clusterNA, fit$'Clinical cluster'))
          names(BL_MR_PD) <- c('cluster1', 'cluster2')
          
          tps <- c("BL", "BR", "PD")
          for (i in 1:length(avail_timepoints)) {
              tp <- tps[avail_timepoints[i]]
              tp_col <- my_data[[as.character(i)]]
              if (!is.null(tp_col)) {
                  BL_MR_PD[[sprintf("VAF_%s", tp)]] <- tp_col
                  BL_MR_PD[[sprintf("is_na_%s", tp)]] <- is.na(tp_col)
              }
          }
          
          total <- BL_MR_PD %>%
              mutate(across(starts_with("VAF"), ~ ifelse(is.na(.), 1e-10, .)))
          total$GENE <- fit$GENE
          total$IS_DRIVER <- fit$IS_DRIVER
          total$IS_ESCAPE <- fit$IS_ESCAPE
          total$significance <- fit$significance
          
          plots <- lapply(1:length(time_combos), function(j) {
              at <- axis_titles[[j]]
              scatter <- ggplot(total, aes(x = .data[[sprintf("VAF_%s", at[1])]], y = .data[[sprintf("VAF_%s", at[2])]], color = factor(significance, levels = my_levels3))) +
                  geom_point(size = 1.5, alpha = 0.5) +
                  scale_color_manual(values = my_colours3) +
                  labs(color = "Clusters", x = at[1], y = at[2]) +
                  theme(
                      panel.background = element_rect(fill = "white", color = "white"),
                      plot.background = element_rect(fill = "white", color = "white"),
                      legend.position = "none",
                      aspect.ratio = 1
                  )
              
              if (sc == "scaled") {
                  scatter <- scatter +
                      scale_x_continuous(breaks = seq(0.1, 1.0, 0.2), labels = signif(mobster:::.log_untransform(seq(0.1, 1.0, 0.2), fit), digits = 2), expand = expansion(mult = c(0.2, 0.01))) +
                      scale_y_continuous(breaks = seq(0.1, 1.0, 0.1), labels = signif(mobster:::.log_untransform(seq(0.1, 1.0, 0.1), fit), digits = 2), expand = expansion(mult = c(0.2, 0.01)))
              }
              
              if (annotate_drivers) {
                  scatter <- scatter +
                      geom_text_repel(data = total %>% filter(IS_DRIVER), aes(x = .data[[sprintf("VAF_%s", at[1])]], y = .data[[sprintf("VAF_%s", at[2])]], label = GENE, fontface = "bold"), inherit.aes = FALSE, nudge_y = 0.025, nudge_x = 0.025, max.overlaps = 1000, size = 2.5, force = 10, segment.color = "grey40", segment.size = 0.1, angle = 0)
              }
              if (annotate_escape) {
                  scatter <- scatter +
                      geom_text_repel(data = total %>% filter(IS_ESCAPE), aes(x = .data[[sprintf("VAF_%s", at[1])]], y = .data[[sprintf("VAF_%s", at[2])]], label = GENE, fontface = "bold"), inherit.aes = FALSE, nudge_y = 0.025, nudge_x = 0.025, max.overlaps = 1000, size = 2.5, force = 10, segment.color = "grey40", segment.size = 0.1, angle = 0)
              }
              return(scatter)
          })
          
          at <- axis_titles[[length(time_combos)]]
          legend <- get_legend(ggplot(total, aes(x = .data[[sprintf("VAF_%s", at[1])]], y = .data[[sprintf("VAF_%s", at[2])]], color = factor(significance, levels = my_levels3))) + geom_point() + scale_color_manual(values = my_colours3) + labs(color = "Clusters"))
          
          # Get the correct table based on plot type
          tables <- mobster:::get_driver_escape_table(fit, cohort)
          if (plot_type=='drivers'){
            table_to_plot = tables$driver_table
          } else {
            table_to_plot = tables$escape_table
          }

          # Combine legend and table into a single plot object
          legend_and_table <- plot_grid(legend, table_to_plot, nrow = 2, rel_heights = c(0.3, 0.7), align = "v")

          # Update the wrapped plots call to include the new object
          wrapped <- wrapped_plots(plots, legend_and_table)
          
          intermed <- plot_grid(subtitle, wrapped$wrapped, ncol = 1, nrow = 2, rel_heights = c(0.025, 1))
          arranged <- plot_grid(title, plot_grid(hist_plot_object, intermed, ncol = 1, nrow = 2, rel_heights = c(0.7, wrapped$nrow)), ncol = 1, nrow = 2, rel_heights = c(0.05, 1))
      }
      
      if (plot_type == 'drivers') {
          pdf_plots_all9[[ploidy_str]] <- arranged
      } else {
          pdf_plots_all10[[ploidy_str]] <- arranged
      }
    }
  }
    
    suptitle <- ggdraw() +
      draw_label(patient_str, fontface = 'bold', x = 0, hjust = 0, size = 40) +
      theme(plot.margin = margin(0, 0, 0, 7))
    
    # Save combined plots for each patient
    plot_height <- 6 * (wrapped$nrow + 1) * length(pdf_plots_all9)
    plot_width <- 15 + (wrapped$ncol) * 5
    
    arranged_all9 <- plot_grid(suptitle, plot_grid(plotlist = pdf_plots_all9, nrow = length(pdf_plots_all9)), nrow = 2, ncol = 1, rel_heights = c(0.05, 1)) + theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
    ggsave(sprintf("%s_scatter_%s_significance_drivers.pdf", patient_str, filtering_str), arranged_all9,
           path = sprintf("%s/clinical_scatters/significance%s%s", output_path, if (sc == "unscaled") { "_untransformed" } else { "" }, if (coverage_filtering) { "_cov_filtered" } else { "" }),
           height = plot_height, width = plot_width, limitsize = FALSE)
    
    arranged_all10 <- plot_grid(suptitle, plot_grid(plotlist = pdf_plots_all10, nrow = length(pdf_plots_all10)), nrow = 2, ncol = 1, rel_heights = c(0.05, 1)) + theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
    ggsave(sprintf("%s_scatter_%s_significance_escape.pdf", patient_str, filtering_str), arranged_all10,
           path = sprintf("%s/clinical_scatters/significance%s%s", output_path, if (sc == "unscaled") { "_untransformed" } else { "" }, if (coverage_filtering) { "_cov_filtered" } else { "" }),
           height = plot_height, width = plot_width, limitsize = FALSE)
  }
}
  
 # pdfunite $(cat file_list.txt) merged.pdf

