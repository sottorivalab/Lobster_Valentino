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
                      no_significance = FALSE,
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
                                                                      ylims = ylims,
                                                                      no_significance = no_significance
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
                            ylims = NULL,
                            no_significance = FALSE
)
{
  MM_scale_a = MM_scale_list$MM_scale_a
  MM_scale_b = MM_scale_list$MM_scale_b
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Plotting data
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Clustering assignments and clusters
  
  x$VAF = if(log_MMscaled==TRUE) x$VAF_ln_scaled else x$VAF_orig 
  if (highlight_NA==FALSE){
    if (no_significance==TRUE){
      plot_data = x %>% dplyr::arrange(cluster)
    } else {
      plot_data = x %>% dplyr::arrange(significance)
    }
  } else {
    plot_data = x %>% dplyr::arrange(clusterNA)
  }
  
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
  if (highlight_NA==FALSE){
    if (no_significance==TRUE){
      hist_pl = ggplot(as_tibble(plot_data),
                       aes(x = VAF[[as.character(t)]],
                           fill = factor(cluster),
                           y = ..count.. / total_count 
                       ))
    } else {
      hist_pl = ggplot(as_tibble(plot_data),
                       aes(x = VAF[[as.character(t)]],
                           fill = factor(significance),
                           y = ..count.. / total_count 
                       ))
    }
  } else {
    hist_pl = ggplot(as_tibble(plot_data),
                     aes(x = VAF[[as.character(t)]],
                         fill = factor(clusterNA),
                         y = ..count.. / total_count 
                     ))
  }
  hist_pl = hist_pl + 
    scale_fill_manual(values = labeled_colours) + 
    geom_histogram(
      # aes(fill)
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
    
    #if ((all(c('GENE','IS_DRIVER') %in% names(as_tibble(plot_data))))&
    #    annotate_drivers){
    #  hist_pl = hist_pl + 
    #    geom_text_repel(
    #      data = as_tibble(plot_data) %>% filter(IS_DRIVER==TRUE),
    #      aes(label = GENE, x = VAF[[as.character(t)]], 
    #          y = rep(0,nrow(as_tibble(plot_data) %>% filter(IS_DRIVER==TRUE))),
    #          fontface = "bold"), 
    #      inherit.aes = FALSE,
    #      max.overlaps = 100,
    #      size = 2.7, # text size
    #      force=5,
    #      nudge_y = 0.01,
    #      force_pull = 1,
    #      segment.color = "grey40", # line colour
    #      segment.size = 0.1,
    #      angle = 90
    #    )  
    #}
  #
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