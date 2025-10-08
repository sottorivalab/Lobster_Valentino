#' Plot a MOBSTER fit.
#'
#' @param x An object of class \code{"dbpmm"}.
#' @param beta_colors A vector of colors that are used to colour the Beta clusters.
#' Colors are used by order for \code{"C1"}, \code{"C2"}, \code{"C3"}, etc. If too
#' few colours are used the \code{rainbow} palette is used instead of \code{beta_colors}.
#' In all cases the colour of the tail is set by the other parameter \code{tail_color}.
#' @param tail_color The colour of the tail cluster, if any.
#' @param cutoff_assignment Parameters passed to run function \code{Clusters} which
#' returns the hard clustering assignments for the histogram plot if one wants to plot
#' only mutations with responsibility above this parameter.
#' @param annotation_extras A dataframe that contains a label column, and a VAF value.
#' The labels will be annotated to the corresponding clusters of the VAF values.
#' @param secondary_axis \code{NULL} to leave the second axis empty, \code{"SSE"} to
#' report the cumulative percentage of SSE error, and \code{"N"} to report the
#' cumulative number of mutations.
#' @param ...
#'
#' @return A ggplot object for the plot.
#'
#' @import sads
#' @import ggplot2
#'
#' @exportS3Method plot dbpmm
#' @export plot.dbpmm
#'
#' @examples
#' data(fit_example)
#' plot(fit_example$best)
#' plot(fit_example$best, secondary_axis = 'SSE')
#' plot(fit_example$best, cutoff_assignment = .7)
#' plot(fit_example$best, beta_colors = c("indianred3", "orange"))
#' plot(fit_example$best, tail_color = 'black')
plot.dbpmm = function(x,
                      cutoff_assignment = 0,
                      beta_colors = RColorBrewer::brewer.pal(n = 9, 'Set1'),
                      tail_color = "gainsboro",
                      annotation_extras = NULL,
                      secondary_axis = NULL,
                      MM_scale_list=NULL, 
                      log_MMscaled=FALSE,
                      transform_axis_back=FALSE,
                      untf = FALSE,
                      ccf = FALSE,
                      highlight_NA = FALSE,
                      annotate_genes = FALSE,
                      xlims = NULL,
                      ylims = NULL,
                      ...)
{
  if(ccf && !transform_axis_back){
    warning('Transforming VAF axis to CCF only available when transform_axis_back == TRUE\n')
  }
  
  mobster:::is_mobster_fit(x)
  
  if (log_MMscaled==TRUE){
    # fit carried out in log space
    
    if (untf==TRUE){ # plot in normal space 
      x$data = x$data %>% dplyr::rename(vaf_scaled = VAF) %>% dplyr::rename(VAF = VAF_orig) 
    }
    
    MM_scale_a = MM_scale_list$MM_scale_a
    MM_scale_b = MM_scale_list$MM_scale_b
    
  } else {
    # fit was not carried out in log space 
    
    if (untf==FALSE){ # display fit in log space
      # transform data into log scale
      MM_scale  = mobster:::get_a_b(log(x$data$VAF))
      MM_scale_a = MM_scale$a 
      MM_scale_b = MM_scale$b 
      
      x$data = x$data %>% mobster:::.log_transform()
      
      # transform scale parameter of pareto
      x$scale = lapply(x$scale, function(x) { MM_scale_a*log(x) + MM_scale_b })
      
      # for use in density functions
      log_MMscaled = TRUE
    } 
  }
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Plotting variables
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  binwidth = 0.01
  domain = seq(0, 1, binwidth)
  
  pi = mobster:::.params_Pi(x)
  
  # Component mixture density and total as \sum_i pi_i f_i
  
  densities = mobster:::template_density(
    x,
    x.axis = domain,
    binwidth = binwidth,
    reduce = TRUE,
    untf = untf, 
    MM_scale_a = MM_scale_a,
    MM_scale_b = MM_scale_b,
    log_MMscaled = log_MMscaled
  )
  
  #if (t==1){
  #  saveRDS(densities,)
  #}
  
  total_densities = lapply(tibble::as_tibble(densities), 
                           function(z){z %>% dplyr::group_by(x) %>%
                               dplyr::summarise(y = sum(y), cluster = 'f(x)')})
  
  #pi = sort(pi)
  #pi = paste0(names(pi), ' ', round(pi * 100, 1), '%', collapse = ', ')
  sor_p = sort(unlist(lapply(x$pi,mean)), decreasing = TRUE)
  sor_p = names(sor_p)
  
  pi = lapply(x$pi[sor_p], function(x){round(x, digits = 2) * 100}) 
  
  hist_pl = lapply(as.list(1:x$n_t),FUN = function(t){individual_plot(t = t, x = x, 
                                                                      domain = domain,
                                                                      densities = densities[[t]],
                                                                      total_densities = total_densities[[t]],
                                                                      pi = as.list(lapply(pi,function(x){x[t]})),
                                                                      binwidth = binwidth,
                                                                      cutoff_assignment = cutoff_assignment,
                                                                      beta_colors = beta_colors,
                                                                      tail_color = tail_color,
                                                                      annotation_extras = annotation_extras,
                                                                      secondary_axis = secondary_axis,
                                                                      MM_scale_list = MM_scale_list, 
                                                                      log_MMscaled = log_MMscaled,
                                                                      transform_axis_back = transform_axis_back,
                                                                      untf = untf,
                                                                      ccf = ccf,
                                                                      highlight_NA = highlight_NA,
                                                                      annotate_genes = annotate_genes,
                                                                      xlims = xlims,
                                                                      ylims = ylims
  )
  })
  
  return(hist_pl)
}

individual_plot <- function(t, x, domain, densities, total_densities, pi, binwidth,
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
                            untf = FALSE,
                            ccf = FALSE,
                            highlight_NA = FALSE,
                            annotate_genes = FALSE,
                            xlims = NULL,
                            ylims = NULL
)
{
  MM_scale_a = MM_scale_list$MM_scale_a
  MM_scale_b = MM_scale_list$MM_scale_b
  
  labels = names(mobster:::.params_Pi(x))
  labels.betas = mobster:::.params_Beta(x)$cluster
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Plotting data
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Clustering assignments and clusters
  #plot_data = mobster:::Clusters(x, cutoff_assignment, log_MMscaled = log_MMscaled,
  #                               MM_scale_a = MM_scale_a, MM_scale_b = MM_scale_b, untf = untf) %>%
  #  dplyr::arrange(cluster) # latent_vars inside Clusters doesn't respect beta cluster permutation
  plot_data = x$data %>% dplyr::arrange(cluster)
  
  if (highlight_NA){ plot_data = posthoc_NA_clusters(plot_data) } else { plot_data$clusterNA = plot_data$cluster }
  
  clusters = sort(unique(plot_data$clusterNA), na.last = TRUE)
  
  # Beta peaks
  if (is.null(MM_scale_b)){
    Beta_peaks = x$Clusters %>%
      dplyr::filter(type == 'Mean', cluster != 'Tail') %>% dplyr::select(-init.value) %>% 
      dplyr::mutate(fit.value = lapply(fit.value,function(x){x[t]}))
  } else { Beta_peaks = get_cluster_mode(x, MM_scale_a, MM_scale_b)}
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Plotting aesthetics
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  # Colours per cluster - tail is reserved
  na_clusters_t =  sapply(possible_NA_clusters_4t(t,x$n_t), function(x)
  {paste('NA_in', paste(x, collapse = "_"), sep = "_")})
  Beta_colors = pio:::nmfy(paste0('C', 1:length(beta_colors)), beta_colors)
  NA_colors = pio:::nmfy(na_clusters_t, na_colors[1:length(na_clusters_t)])
  
  #if ('gene' %in% names(as_tibble(plot_data))){
  #  Driver_colors = pio:::nmfy(unique(plot_data$gene), driver_colors[1:length(unique(plot_data$gene))])
  #  MOBSTER_CLUSTER_COLORS = c(Beta_colors, NA_colors, `Tail` = tail_color, Driver_colors)
  #} else {
  MOBSTER_CLUSTER_COLORS = c(Beta_colors, NA_colors, `Tail` = tail_color)
  #}
  
  if (all(clusters %in% names(MOBSTER_CLUSTER_COLORS))|(highlight_NA))
  {
    MOBSTER_CLUSTER_COLORS = MOBSTER_CLUSTER_COLORS[names(MOBSTER_CLUSTER_COLORS) %in% clusters]
  }
  else
  {
    available =  clusters[clusters %in% names(MOBSTER_CLUSTER_COLORS)]
    missing =  clusters[!(clusters %in% names(MOBSTER_CLUSTER_COLORS))]
    
    warning(
      "You did not pass enough input colours, adding a gray colour\n",
      'Available: ', paste0(available, collapse = ', '), '\n',
      'Missing: ', paste0(missing, collapse = ', ')
    )
    
    MOBSTER_CLUSTER_COLORS = c(MOBSTER_CLUSTER_COLORS,
                               pio:::nmfy(paste(missing), rep("gray", length(missing)))
    )
    # MOBSTER_CLUSTER_COLORS = rainbow(length(clusters))
  }
  
  # Plot title
  plot_title = ifelse(all(is.null(x$description)), 'lobster fit', x$description)
  
  # Plot substitle
  plot_subtitle = NULL
  pi_label = paste0(pi, '% [', names(pi), ']', collapse = ', ')
  
  #plot_subtitle = paste0("n = ", x$N, "; ", pi_label)
  plot_subtitle = paste0("n = ", nrow(plot_data[!is.na(plot_data$VAF[[as.character(t)]]), ]), "; ", pi_label)
  
  # Plot caption
  plot_caption = NULL
  
  conv.steps = length(x$all.NLL)
  conv.epsilon = 0
  if (conv.steps >= 2){
    #conv.epsilon = abs(rev(x$all.NLL)[1] - rev(x$all.NLL)[2])
    # changed to give fraction relative to total decrease between first and last
    conv.epsilon = abs((rev(x$all.NLL)[1] - rev(x$all.NLL)[2])/(rev(x$all.NLL)[1]-x$all.NLL[1]))
  }
  conv.epsilon = formatC(conv.epsilon, format = "e", digits = 0)
  
  plot_caption = bquote(
    .(x$fit.type) *
      " (" * .(conv.steps) ~ 'steps; fractional ' * epsilon ~ '=' ~ .(conv.epsilon) *
      '), ' * z["nk"] * ' > ' * .(cutoff_assignment)
  )
  
  total_count <- nrow(plot_data[complete.cases(plot_data$VAF[[as.character(t)]]), ])
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Main ggplot object is the histogram
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  hist_pl = ggplot(as_tibble(plot_data),
                   aes(x = VAF[[as.character(t)]],
                       fill = factor(clusterNA),
                       y = ..count.. / total_count #sum(..count..)
                   )) +
    geom_histogram(
      # aes(fill)
      alpha = .6,
      color = NA,
      position = 'stack',#'identity',
      binwidth = binwidth
    ) +
    scale_fill_manual(values = MOBSTER_CLUSTER_COLORS) +
    geom_vline(
      xintercept = min(x$data$VAF,na.rm = TRUE),
      colour = 'black',
      linetype = "longdash",
      size = .3
    ) +guides(fill = guide_legend(title = "Cluster")) + 
    geom_line(data = densities, 
              aes(
                y = y,
                x = x,
                color = factor(cluster)
              ),
              size = 1, inherit.aes = FALSE) +
    geom_line(
      data = total_densities,# %>% mutate(value = value + max(.$value, na.rm = TRUE) * 0.02), # changed
      aes(y = y, x = x),
      color = 'black',
      alpha = .8,
      size = .5,
      linetype = 'dashed',
      inherit.aes = FALSE
    ) +
   scale_color_manual(values = MOBSTER_CLUSTER_COLORS) +
    guides(color = FALSE) +
    labs(
      title = bquote(.(plot_title)),
      subtitle = plot_subtitle,
      caption = plot_caption,
      x = ifelse((!ccf),"Observed VAF","Observed CCF"),
      y = "Density"
    ) +
    mobster:::my_ggplot_theme() +
    coord_cartesian(clip = 'off')
   # theme( # error here
   #   panel.background = element_rect(fill = 'white'),
   #   plot.caption = element_text(color = ifelse(x$status, "darkgreen",  "red"))
   # ) #+
  
  #if ('gene' %in% names(as_tibble(plot_data))){
  if ((all(c('GENE','IS_DRIVER') %in% names(as_tibble(plot_data))))&
      annotate_genes){
    hist_pl = hist_pl + 
      geom_text_repel(
        data = as_tibble(plot_data) %>% filter(IS_DRIVER==TRUE),
        aes(label = GENE, x = VAF[[as.character(t)]], 
            y = rep(0,nrow(as_tibble(plot_data) %>% filter(IS_DRIVER==TRUE))),
            fontface = "bold"), 
        inherit.aes = FALSE,
        nudge_y = -0.1,
        max.overlaps = 100,
        size = 2, # text size
        force=1,
        segment.color = "grey20", # line colour
        segment.size = 0.1,
        angle = 90
      )  
    # ) 
  }
  # if (is.null(MM_scale_b)){
  # hist_pl = hist_pl +  geom_vline(data = Beta_peaks,
  #                      aes(xintercept = unlist(fit.value),
  #                          color = factor(cluster)),
  #                      linetype = "longdash")
  # } else {
  #   peaks_df <- data.frame(cluster = names(Beta_peaks), value = unlist(Beta_peaks))  
  #   hist_pl = hist_pl +  geom_vline(data = peaks_df,
  #                        aes(xintercept = value,
  #                            color = factor(cluster)),
  #                        linetype = "longdash")
  # }
  
  #if (!is.null(xlims)){
  #  print(here)
  #  hist_pl = hist_pl + xlim(xlims)
  #}
  
  if (transform_axis_back){
    my_seq = seq(0.1,1,0.15)
    if (ccf){
      hist_pl = hist_pl +
        scale_x_continuous(breaks = my_seq,  
                           labels = signif(0.5*ccf_to_vaf(.log_untransform(my_seq,x$data),x,t), digits=2)) 
    } else {
      hist_pl = hist_pl +
        scale_x_continuous(breaks = my_seq,  
                           labels = signif(mobster:::.log_untransform(my_seq,x$data), digits=2)) 
    }
  }  
  if (untf){
   hist_pl = hist_pl + xlim(0, 0.5) 
  }
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Dashed SSE behind the overall density (optional), scaled to percentage
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # LEAVE TILL LATER
  #if (!is.null(secondary_axis))
  #{
  #  if (secondary_axis == "SSE")
  #  {
  #    error = suppressWarnings(mobster:::.compute_fit_sqerr(x, binning = binwidth,
  #                                                          log_MMscaled = log_MMscaled, 
  #                                                          MM_scale_a = MM_scale_a,
  #                                                          MM_scale_b = MM_scale_b))
  #    
  #    m = max(ggplot_build(hist_pl)$data[[1]]$y)
  #    
  #    me = max(error$cum.y)
  #    error = error %>% dplyr::mutate(cum.y = (cum.y / me) * m)
  #    
  #    hist_pl = hist_pl +
  #      geom_line(
  #        data = error,
  #        aes(y = cum.y, x = x),
  #        color = 'darkgray',
  #        alpha = 1,
  #        size = .2,
  #        linetype = 'dashed',
  #        inherit.aes = FALSE
  #      ) +
  #      scale_y_continuous(sec.axis = sec_axis(~ . / m * 100, name = "SSE [cumulative %]"))
  #  }
  
  # if(secondary_axis == "N")
  # {
  #   xbins = hist(x$data$VAF, breaks = domain, plot = F)$counts
  #   xcumc = cumsum(xbins)
  #
  #   brk = seq(0, x$N, by = x$N/4)
  #
  #   # hist_pl =
  #     hist_pl +
  #     geom_line(
  #       data =  data.frame(x = domain[-1], y = (xcumc/x$N) * m),
  #       aes(y = y, x = x),
  #       color = 'darkgray',
  #       alpha = 1,
  #       size = .2,
  #       linetype = 'dashed',
  #       inherit.aes = FALSE
  #     ) +
  #     scale_y_continuous(
  #       sec.axis = sec_axis( ~ . ,
  #                            name = "Counts (n)",
  #                            breaks = brk / x$N * m,
  #                            labels = brk
  #                            )
  #       )
  # }
  
  # }
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Extra annotations (param) plus drivers
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # hist_pl = suppressWarnings(
  #  mobster:::add_extra_plot_annotations(x, annotation_extras, base_plot = hist_pl)
  #)
  
  # setting aspect ratio
  #hist_pl = hist_pl + coord_fixed(ratio = ifelse(annotate_genes,1,0.75))
  
  # Get the space needed for labels
  current_ylims <- ggplot_build(hist_pl)$layout$panel_params[[1]]$y.range
  
  # Automatically adjust ylim
  if (annotate_genes & is.null(ylims)){
    hist_pl = hist_pl + ylim(c(current_ylims[1] - 0.15, current_ylims[2]))
  } else if (!is.null(ylims)){
    hist_pl = hist_pl + ylim(ylims)
  }
  
  return(hist_pl)
}