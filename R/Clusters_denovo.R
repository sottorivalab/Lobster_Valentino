#' Assign new observations to the clusters inside a MOBSTER fit.
#' 
#' @description For a new set of observations in the MOBSTER input format,
#' the parameters and clusters of a MOBSTER fit are used to determine the
#' assignments of the new points. MOBSTER density function is used to determine
#' which cluster has the highest density for each observation, and the cluster
#' label is computed accordingly. 
#'
#' @param x A MOBSTER fit.
#' @param y A MOBSTER input dataset, which has to have a VAF numeric column without NAs.
#'
#' @return The data in \code{y} is augmented with a colum per mixture component reporting
#' the corresponding density value. A final colum `cluster` is also added reporting the
#' component name for the hard clustering assignment of the point.
#' 
#' @export
#'
#' @examples
#' library(ggplot2)
#' data('fit_example', package = 'mobster')
#' 
#' # Generate some randome numbers and assign them to the most likely mixture component
#' new_assignments = Clusters_denovo(fit_example$best, data.frame(VAF = runif(1000))) 
#' print(new_assignments)
#' 
#' # Plot a histogram coloured according to the clusters
#' ggplot(new_assignments, aes(VAF, fill = cluster)) + geom_histogram(binwidth = 0.01)

Clusters_denovo = function(x, y)
{
  if (nrow(y)>0){
    mobster:::is_mobster_fit(x)
    #mobster:::is_mobster_input_dataset(y) # update later 
    
    z_nkt     = array(0, c(nrow(y),x$K,x$n_t))  
    pdf.w     = array(0, c(nrow(y),x$K,x$n_t))  
    
    for (k in 1:x$K){
      for (t in 1:x$n_t){
        pdf.w[,k,t] = mobster::ddbpmm(x,
                                      data = y$VAF,
                                      components = k,
                                      a = x$a,
                                      b = x$b,
                                      pi = x$pi,
                                      shape = x$shape,
                                      scale = x$scale,
                                      log = TRUE,
                                      sigma = x$sigma[t],
                                      time_point = t,
                                      untf = TRUE
        )
      }
    }
    
    # for calculating pi (for each t)
    Z       = apply(pdf.w, c(1,3), mobster:::.log_sum_exp)
    Z   = array(Z, dim = c(dim(z_nkt)[1],dim(z_nkt)[3])) # protect against dropping one-element dims 
    
    z_nkt   = sweep(pdf.w,MARGIN=c(1,3),Z,FUN="-")  # x$pdf.w - Z 
    z_nkt_exp   = apply(z_nkt, c(2,3), exp) # Exponentiate to get actual probabilities
    z_nkt   = array(z_nkt_exp, dim = c(dim(z_nkt))) # protect against dropping one-element dims 
    
    resp_prods = apply(z_nkt, c(1,2), function(z){mean(z,na.rm = TRUE)}) # multiply responsibilities across timepoints (CHECK!)
    y$cluster_denovo =  unlist(apply(resp_prods, 1,
                                     function(z) {names(x$pi)[which(z == max(z, na.rm = TRUE))[1]]}))
    
  }
      
  return(y)
}
