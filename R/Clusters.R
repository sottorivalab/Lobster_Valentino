#' Return the data with the hard clustering assigments.
#' 
#' @description Extract the clustering assignments, and return the data.
#' Assignments can be computed subsetting the mutations that have posterior
#' values (reponsibilities) below a certain cutoff (default \code{0} - all 
#' assigments); non-assigned mutations have \code{NA} as cluster label.
#'
#' @param x A MOBSTER fit.
#' @param cutoff_assignment Cutoff to compute hard clustering assignments.
#'
#' @return The data stored in \code{x$data} with a column \code{label} reporting
#' the assigned cluster, or \code{NA} if the maximum cluster probability is below
#' the threshold value \code{cutoff_assignment}.
#' 
#' @export
#'
#' @examples
#' data('fit_example', package = 'mobster')
#' Clusters(fit_example$best)
#' 
#' # Add some cutoff to filter assignments
#' Clusters(fit_example$best, cutoff_assignment = .8)
Clusters = function(x, cutoff_assignment = 0, log_MMscaled = FALSE, MM_scale_a = NULL, MM_scale_b = NULL, untf = FALSE)
{
  browser()
  is_mobster_fit(x)
  
  lv = mobster:::latent_vars(x, log_MMscaled = log_MMscaled, MM_scale_a = MM_scale_a, MM_scale_b = MM_scale_b, untf = untf)

  # Maximum density form LV
  #maxima = apply(lv$z_nkt, 1, max)
  
  # LV and posterior hard clusters
  assignments = mobster:::latent_vars_hard_assignments(lv)
  
  #assignments[maxima < cutoff_assignment] = NA
  
  x$data$cluster = assignments
  
  dplyr::bind_cols(x$data, lv$z_nkt %>% as_tibble())
}