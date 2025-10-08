#' Density function for MOBSTER.
#'
#' @description The density function for a MOBSTER mixture. The function
#' allows to compute the density from the data and parameters stored
#' inside an object of class \code{dbpmm}, or from
#' a custom set of parameters and data that can be passed as input. The
#' function can also compute the density only for a subset of components
#' of the mixture.
#'
#' @param x An object of class \code{"dbpmm"}.
#' @param data Data to compute the (log)-likelihood, if this is null the data stored
#' inside \code{x$data} is used.
#' @param components A vector specifying which components should be used to compute
#' the density. Position 1 is reserved for the tail (regardless that is fit or not
#' to the data), the other positions refer to Beta components.
#' @param log Boolean value to select the log-likelihood, or the likelihood.
#' @param a Vector of parameters for the Beta components; if this is null values stored
#' inside \code{x} are used.
#' @param b Vector of parameters for the Beta components; if this is null values stored
#' inside \code{x} are used.
#' @param pi Mixing proportions for the mixture; if this is null values stored
#' inside \code{x} are used.
#' @param scale Scale of the power law; if this is null values stored
#' inside \code{x} are used.
#' @param shape Shape of the power law; if this is null values stored
#' inside \code{x} are used.
#'
#' @return The density for the data. Notice that to compute the negative log-likelihood
#' used during fit one needs to use \code{log = TRUE} and change sign.
#' 
#' @export
#' 
#' @importFrom  sads dpareto
#'
#' @examples 
#' library(ggplot2)
#' data('fit_example', package = 'mobster')
#' 
#' # Use the full mixture, and its internal data
#' ddbpmm(fit_example$best)
#' 
#' # Use only some of the mixture components, and pass some data
#' ddbpmm(fit_example$best, data = .4, components = 1)
#' 
#' # An internal function to get f(x) with x the [0,1] range.
#' ggplot(mobster:::template_density(fit_example$best, reduce = TRUE),
#'        aes(x = x, y = y, color = cluster)) + geom_line()

erf <- function(x) {2*pnorm(x*sqrt(2))-1} # taken from SO

gauss_norm_func = function(x, sigma, MMa, MMb){
  gn = MMa*sigma*sqrt(2*pi)*(1 - erf(-exp((x-MMb)/MMa)/(sigma*sqrt(2))))/2
  return(gn)
}

gauss_norm_untf = function(x, sigma){
  # for convolution with x between 0 and 1
  gn = (erf((1-x)/(sigma*sqrt(2))) - erf(- x /(sigma*sqrt(2))))/2
  return(gn)
}

gauss_MMscaled_wderiv = function(y,x,sigma,MMa,MMb){
  gf = exp((y-MMb)/MMa)*exp(-0.5*((exp((y-MMb)/MMa)-exp((x-MMb)/MMa))/sigma)**2) 
  return(gf)
}

gauss_MMscaled_wderiv_normed = function(y,x,sigma,MMa,MMb){
  return(exp(log(gauss_MMscaled_wderiv(y,x,sigma,MMa,MMb)) - log(gauss_norm_func(x,sigma,MMa,MMb))))
}

gauss_untf = function(y,x,sigma){
  return((1/(sigma*sqrt(2*pi)))*exp(-0.5*((y-x)/sigma)**2))
}

d_expMMtail = function(data, shape, scale, MMa){
  xstar = scale
  multiplier = crone::heaviside(data, x0 = xstar) + .Machine$double.xmin
  density = exp(log(multiplier)+log(shape) + (shape * (scale-data) / MMa) - log(MMa)) # checked
  return(density)
}

# unnormalised Exponentiated MM-scaled beta density
d_expMMbeta_unnorm = function(data, shape1, shape2, MMa, MMb){
  
  density_unnorm = (1/MMa)*exp((shape1/MMa)*(data-MMb) + (shape2-1)*log((1 - exp((1/MMa)*(data-MMb)))))
  
  return(density_unnorm)
}

d_expMMtail_convolved = function(data, shape, scale, MMa, MMb, sigma, log = FALSE){
  # here, scale should be transformed value
  integrand <- function(x,y) {exp(log(gauss_MMscaled_wderiv_normed(y,x,sigma,MMa,MMb))+log(d_expMMtail(x,shape,scale,MMa)))}
  integrator <- function(y) {
    return(integrate(integrand,scale,MMb,abs.tol = sqrt(.Machine$double.eps),y=y,stop.on.error=FALSE)$value)} 
  
  integral <- Vectorize(integrator)
  #dat_integral <- parallel::pvec(data, integral, mc.cores=7) 
  dat_integral = integral(data)
  
  if (log==TRUE){
    tryCatch({dat_integral = log(dat_integral)},
             error = function(e){dat_integral = rep(NA,length(dat_integral))})
  }
  return(dat_integral)
}

d_expMMbeta_convolved = function(data, shape1, shape2, MMa, MMb, sigma, log = FALSE){
  
  beta_MMscaled_partition = beta(shape1,shape2)
  integrand <- function(x,y) {exp( log(gauss_MMscaled_wderiv_normed(y,x,sigma,MMa,MMb))+
                                     log(d_expMMbeta_unnorm(x,shape1,shape2,MMa,MMb)) - log(beta_MMscaled_partition))}
  integrator <- function(y) {
    return(integrate(integrand,-Inf,MMb,abs.tol = sqrt(.Machine$double.eps),y=y,stop.on.error=FALSE)$value)} 
  
  integral <- Vectorize(integrator)
  dat_integral = integral(data)
  #dat_integral <- parallel::pvec(data, integral, mc.cores=7)
  
  if (log==TRUE){
    tryCatch({dat_integral = log(dat_integral)},
             error = function(e){dat_integral = rep(NA,length(dat_integral))})
  }
  
  return(dat_integral)
}

d_pareto_convolved <- function(data, shape, scale, sigma, log = FALSE){
  integrand <- function(x,y) {exp( log(gauss_untf(y,x,sigma)) +
                                     sads::dpareto(x,shape = shape,scale = scale,log = TRUE) -
                                     log(gauss_norm_untf(x,sigma))) }
  integrator <- function(y) {
    return(integrate(integrand,0.0,1,abs.tol = sqrt(.Machine$double.eps),y=y,stop.on.error=FALSE)$value)} 
  
  integral <- Vectorize(integrator)
  #dat_integral <- parallel::pvec(data, integral, mc.cores=7) 
  dat_integral = integral(data)
  
  if (log==TRUE){dat_integral = log(dat_integral)}
  return(dat_integral)
}

d_beta_convolved <- function(data, a, b, sigma, log = FALSE){
  integrand <- function(x,y) { exp( log(gauss_untf(y,x,sigma)) +
                                      dbeta(x,shape1 = a,shape2 = b,log = TRUE) -
                                      log(gauss_norm_untf(x,sigma)) ) }
  integrator <- function(y) {
    return(integrate(integrand,0,1,abs.tol = sqrt(.Machine$double.eps),y=y,stop.on.error=FALSE)$value)} 
  
  integral <- Vectorize(integrator)
  #dat_integral <- parallel::pvec(data, integral, mc.cores=7) 
  dat_integral = integral(data)
  
  if (log==TRUE){dat_integral = log(dat_integral)}
  return(dat_integral)
}

density_NA_handler = function(data, my_func, ...){
  # handles NA in data while preserving vectorisation
  NA_mask = is.na(data)
  nonNA_idx = which(!NA_mask)
  term = rep(NA,length(data))
  if (length(nonNA_idx)>0){
    filtered_data = data[!NA_mask]
    term[nonNA_idx] = my_func(filtered_data,...) 
  }
  
  return(term)
}

ddbpmm = function(x,
                  data = NULL,
                  components = 1:x$K,
                  a = NULL,
                  b = NULL,
                  pi = NULL,
                  scale = NULL,
                  shape = NULL,
                  log = TRUE,
                  log_MMscaled = FALSE,
                  MM_scale_a = NULL, 
                  MM_scale_b = NULL, 
                  sigma = NULL,
                  untf = FALSE,
                  time_point = NULL, 
                  time = NULL) 
{
  # implementation for user-defined 'time' parameter incomplete - 
  is_mobster_fit(x)
  stopifnot(length(components) > 0)
  
  # Get parameters and data if these are not passed explicitly
  if (all(is.null(data)))
    data = x$data$VAF
  
  if (is.null(scale))
    scale = unlist(.params_Pareto(x)$Scale)
  if (is.null(shape))
    shape = unlist(.params_Pareto(x)$Shape)
  else{ shape = unlist(shape) }
  
  if (all(is.null(a))){
    if (x$K>1){
      a = .params_Beta(x)$a
    }
  }
  if (all(is.null(b))){
    if (x$K>1){
      b = .params_Beta(x)$b
    }
  }
  
  if (all(is.null(pi)))
    pi = .params_Pi(x)
  
  if (is.null(sigma))
    sigma = x$sigma
  
  n_t = x$n_t
  
  scale = scale %>% unlist()
  if (untf & log_MMscaled){scale = mobster:::.log_untransform(scale,x$data)}
  
  # Mask what we do not compute
  mask = (1:x$K) %in% components
  # time mask
  it = 1:n_t
  ts = x$data$times[1,]
  
  if (!is.null(time_point)){ #time point provided
    tmask = rep(0,n_t)
    tmask[time_point] = 1
  } 
  else if (!is.null(time)){ # time in units provided
    tmask = c(1)
    it = 1
    ts = list(time) 
  } 
  else {tmask = rep(1,n_t)} # if no specific time or time-point provided, sum logLik over all timepoints
  
  log_pi = array(sapply(pi,log), dim = c(n_t, x$K))
  
  logLik = 0
  for (t in it){
    if (tmask[t]==1){
      # Tail
      if (mask[1])
      {
        if (x$fit.tail){
          if (untf){
            logLik = log_pi[t,1] + density_NA_handler(data[[t]],d_pareto_convolved,
                                                      shape = shape[t],
                                                      scale = scale,
                                                      sigma = sigma,
                                                      log = TRUE) * mask[1]
          } else if (!log_MMscaled){ # normal Pareto tail 
            logLik = log_pi[t,1] + sads::dpareto(data[[t]],
                                                 shape = shape[t],
                                                 scale = scale,
                                                 log = TRUE) * mask[1]
          } else { # exponential (min-max scaled) tail
            if (is.null(MM_scale_a))
            {stop("MM_scale_a must be provided if log_MMscaled is set to TRUE")}
            logLik = log_pi[t,1] + density_NA_handler(data[[t]],d_expMMtail_convolved,
                                                      shape = shape[t],
                                                      scale = scale,
                                                      MMa = MM_scale_a,
                                                      MMb = MM_scale_b,
                                                      sigma = sigma,
                                                      log = TRUE
            ) * mask[1]
          } }
        else { # Otherwise -Inf is the min achievable
          if (all(class(data)=="numeric")) # this construction is horrible
            logLik = rep(-Inf, length(data)) 
          else
            logLik = rep(-Inf, nrow(data)) 
        }
      }
      
      # Beta components
      if (x$K>1){
        log_pi_t = log_pi[t,2:x$K]
        
        for (k in 1:length(log_pi_t))
          if (mask[k + 1]){
            if (untf){
              logLik = logLik + log_pi_t[k] + density_NA_handler(data[[t]],d_beta_convolved,
                                                                 a = a[[k]][t],
                                                                 b = b[[k]][t],
                                                                 sigma = sigma,
                                                                 log = TRUE)
            } else if (is.null(MM_scale_b)){
              logLik = logLik + log_pi_t[k] + dbeta(data[[t]],
                                                    shape1 = a[[k]][t],
                                                    shape2 = b[[k]][t],
                                                    log = TRUE)
            } else {
              logLik = logLik +
                log_pi_t[k] + density_NA_handler(data[[t]],d_expMMbeta_convolved,
                                                 shape1 = a[[k]][t],
                                                 shape2 = b[[k]][t],
                                                 MMa = MM_scale_a,
                                                 MMb = MM_scale_b,
                                                 sigma = sigma,
                                                 log = TRUE)
            }
          }
      }
    }
  }
  
  if (!log)
    logLik = exp(logLik)
  
  return(logLik)
}


template_density = function(x,
                            x.axis = seq(0, 1, 0.01),
                            init = FALSE,
                            binwidth = 0.01,
                            reduce = FALSE,
                            log_MMscaled = FALSE,
                            MM_scale_a = NULL,
                            MM_scale_b = NULL,
                            set_a = NULL,
                            set_b = NULL,
                            untf = FALSE,
                            time_point = NULL,
                            time = NULL)
{
  if (!is.null(time)){
    stop('time parameter not implemented yet')
  }
  
  if (!is.null(time_point)){
    t = time_point
  } else {t = 1:x$n_t}
  
  data = as_tibble(replicate(x$n_t, x.axis))
  
  labels = names(mobster:::.params_Pi(x))
  values2 = list()
  for (tp in t){
    values = lapply(1:x$K,
                    # Component wise
                    FUN = function(w){
                      data.frame(
                        cluster = labels[w],
                        x = x.axis,
                        y = mobster::ddbpmm(
                          x,
                          data = data,
                          components = w,
                          # init = init,
                          log = FALSE,
                          log_MMscaled=log_MMscaled,
                          MM_scale_a=MM_scale_a,
                          MM_scale_b=MM_scale_b,
                          a = set_a, # for debugging/visualisation purposes
                          b = set_b,
                          untf = untf,
                          time_point = tp,
                          time = time,
                          sigma = x$sigma[tp]
                        ))}) 
    
    names(values) = labels
    
    # remove tail if not required, which is entry #1
    if (!x$fit.tail)
      values = values[2:x$K]
    
    # Scale density wrt binwidth
    values = lapply(values, function(w) {
      w$y = w$y * binwidth
      w
    })
    
    # Reduce if required
    if (reduce)
      values = Reduce(rbind, values)
    
    values2 = append(values2,list(values))
  }
  
  names(values2) = sapply(1:x$n_t,function(x) sprintf('t_%d',x))
  
  values2
}



# UNUSED FUNCTION: numerical nlogLik for a Pareto component (the actual solution is analytical)
# a, b parameters of the mixture, k the component ID, and x data
# Creates a functional that depends only on a and b
# .NLLParetoMix = function(x, z_nk, pi, scale)
# {
#   f = function(shape) # NLL
#   {
#     # z[x,k] * { log pi[k] + log Pareto(x | shape, scale) }
#     -1 * z_nk[, 1] %*% (log(pi[1]) + sads::dpareto(x, shape = shape, scale = scale, log = TRUE))
#   }
#
#   return(f)
# }


#' Generate a random sample from a MOBSTER model.
#' 
#' @description This function can be used to generate data from a MOBSTER
#' model, or from a custom mixture. The principle is the same of function
#' \code{\link{ddbpmm}}.
#'
#' @param x An object of class \code{"dbpmm"}.
#' @param n Number of samples required (i.e., points).
#' @param a Vector of parameters for the Beta components; if this is null values stored
#' inside \code{x} are used.
#' @param b Vector of parameters for the Beta components; if this is null values stored
#' inside \code{x} are used.
#' @param pi Mixing proportions for the mixture; if this is null values stored
#' inside \code{x} are used.
#' @param scale Scale of the power law; if this is null values stored
#' inside \code{x} are used.
#' @param shape Shape of the power law; if this is null values stored
#' inside \code{x} are used.
#' @param tail.cutoff Because the Pareto Type I power law has support
#' over all the positive real line, tail values above \code{tail.cutoff}
#' are rejected (and re-sampled). 
#'
#' @return n samples from the mixture
#' 
#' 
#' @export
#'
#' @examples
#' library(ggplot2)
#' # 1 Beta component at 0.5 mean (symmetrical Beta)
#' a = b = 50
#' names(a) = names(b) = "C1"
#'
#' # 60% tail mutations
#' pi = c('Tail' = .6, 'C1' = .4)
#'
#' # Sample
#' v = data.frame(x = rdbpmm(
#'    x = NULL,
#'    n = 1000,
#'    a = a,
#'    b = b,
#'    pi = pi,
#'    shape = 2,
#'    scale = 0.05
#'  ))
#' ggplot(v, aes(x)) + geom_histogram(binwidth = 0.01)
#' 
#' # Or use the parameters of a model available
#' data('fit_example', package = 'mobster')
#' 
#' v = data.frame(x = rdbpmm(x = fit_example$best, n = 1000))
#' ggplot(v, aes(x)) + geom_histogram(binwidth = 0.01)
rdbpmm = function(x,
                  a = NULL,
                  b = NULL,
                  pi = NULL,
                  shape = NULL,
                  scale = NULL,
                  n = 1,
                  tail.cutoff = 1)
{
  
  # Get parameters if these are not passed explicitly
  if (all(is.null(scale)) | all(is.null(shape)))
  {
    scale = .params_Pareto(x)$Scale
    shape = .params_Pareto(x)$Shape
  }
  
  if (all(is.null(a)) | all(is.null(b)))
  {
    Betas = .params_Beta(x)
    a = Betas$a
    b = Betas$b
    names(a) = names(b) = Betas$cluster
  }
  
  if (all(is.null(pi)))
    pi = .params_Pi(x)
  
  ############################ Sample mixing prop.
  pi.samples = sample(names(pi),
                      size = n,
                      prob = pi,
                      replace = TRUE)
  
  ############################ Tail samples
  n.tail = table(pi.samples)['Tail']
  sample.tails = NULL
  
  if (!is.na(n.tail) & n.tail > 0)
  {
    N = n.tail
    
    # as the model is an approximation to the VAF we reject values > 1.
    repeat {
      # Get `N` Pareto Type I samples
      new.sample.tails = sads::rpareto(N,
                                       shape = shape,
                                       scale = scale)
      
      # Store only values <= 1, which might be less than N
      new.sample.tails = new.sample.tails[new.sample.tails <= tail.cutoff]
      
      sample.tails = c(sample.tails, new.sample.tails)
      
      # Stop if we accumulated `n.tail` samples
      if (length(sample.tails) == n.tail)
        break
      
      
      # Otherwise repeat drawing N-length(sample.tails)
      N = N - length(new.sample.tails)
    }
  }
  
  ############################ Beta samples
  n.Betas = table(pi.samples)[names(table(pi.samples)) != 'Tail']
  
  sample.Betas = sapply(names(n.Betas),
                        function(w)
                        {
                          rbeta(n.Betas[w], shape1 = a[w], shape2 = b[w])
                        })
  
  return(c(unlist(sample.Betas), unlist(sample.tails)))
}

# USED FOR NUMERICAL SOLUTION OF THE MLE ESTIMATORES FOR EACH BETA COMPONENT
# a, b parameters of the mixture, k the component ID, and x data
# Creates a functional that depends only on a and b
#.NLLBetaMix = function(k, x, z_nk, pi)
#{
#  f = function(a,b) # NLL
#  { return(-1 * z_nk[, k] %*% (log(pi[k]) + dbeta(x, shape1 = a, shape2 = b, log = TRUE)))}
#  return(f)
#}

.NLLBetaMixL = function(k, x, z_nkt, pi, t)
{
  f = function(a,b) # NLL
  { 
    sum = - (z_nkt[,k,t][!is.na(z_nkt[,k,t])] %*% (log(pi[[k]][t]) + dbeta(x[[t]][!is.na(z_nkt[,k,t])], shape1 = a, shape2 = b, log = TRUE)))
    return(sum)
  }
  return(f)
}

.NLLBetaMixG = function(k, x, z_nkt, pi, t)
{
  f = function(y) # NLL
  { 
    sum = - (z_nkt[,k,t][!is.na(z_nkt[,k,t])] %*% (log(pi[[k]][t]) + dbeta(x[[t]][!is.na(z_nkt[,k,t])], shape1 = y[1], shape2 = y[2], log = TRUE)))
    return(sum)
  }
  return(f)
}

.NLLBetaMixL_convolved = function(k, x, z_nkt, pi, sigma, t)
{
  f = function(a,b) # NLL
  { 
    sum = - (z_nkt[,k,t][!is.na(z_nkt[,k,t])] %*% (log(pi[[k]][t]) + d_beta_convolved(x[[t]][!is.na(z_nkt[,k,t])], a = a, b = b, sigma = sigma[t], log = TRUE)))
    return(sum)
  }
  return(f)
}

.NLLBetaMixG_convolved = function(k, x, z_nkt, pi, sigma, t)
{
  f = function(y) # NLL
  { 
    sum = - (z_nkt[,k,t][!is.na(z_nkt[,k,t])] %*% (log(pi[[k]][t]) + d_beta_convolved(x[[t]][!is.na(z_nkt[,k,t])], a = y[1], b = y[2], sigma = sigma[t], log = TRUE)))
    return(sum)
  }
  return(f)
}

# without fitting the scale
.NLL_Pareto = function(x, z_nkt, pi, scale, t)
{
  f = function(shape){ # NLL
    sum = - (z_nkt[,1,t][!is.na(z_nkt[,1,t])]) %*% (log(pi[[1]][t]) + sads::dpareto(x[[t]][!is.na(z_nkt[,1,t])], shape = shape, scale = scale, log = TRUE)) 
    return(sum)
  }
  return(f)
}

# without fitting the scale
.NLL_Pareto_convolved = function(x, z_nkt, pi, scale, sigma, t)
{
  f = function(shape){ # NLL
    sum = - (z_nkt[,1,t][!is.na(z_nkt[,1,t])]) %*% (log(pi[[1]][t]) + d_pareto_convolved(x[[t]][!is.na(z_nkt[,1,t])], 
                                                                                         shape = shape, scale = scale, sigma = sigma[t], log = TRUE))
    
    return(sum)
  }
  return(f)
}

# 'transformed' densities are used for visualisation/plots
.NLL_ExpMMBetaMixL_convolved = function(k, x, z_nkt, pi, MMa, MMb, sigma, t) # pass x[[t]] here
{
  
  f = function(a,b) # NLL
  { 
    if ((a<0)|(b<0)|(a>75)|(b>350)){ return(Inf) } 
    
    tryCatch({
      sum = - (z_nkt[,k,t][!is.na(z_nkt[,k,t])] %*% (log(pi[[k]][t]) + 
                                                       d_expMMbeta_convolved(x[!is.na(z_nkt[,k,t])], shape1 = a, shape2 = b,
                                                                             MMa = MMa, MMb = MMb, sigma=sigma[t], log = TRUE)))
    },
    error = function(e){sum = Inf})
    return(sum)
  }
  return(f)  
}

.NLL_ExpMMBetaMixG_convolved = function(k, x, z_nkt, pi, MMa, MMb, sigma, t) # pass x[[t]] here
{
  
  f = function(y) # NLL
  { 
    a = y[1]
    b = y[2]  
    
    if ((a<0)|(b<0)|(a>75)|(b>350)){ return(Inf) } 
    
    sum = 0
    tryCatch({
      sum = sum - (z_nkt[,k,t][!is.na(z_nkt[,k,t])] %*% (log(pi[[k]][t]) + 
                                                           d_expMMbeta_convolved(x[!is.na(z_nkt[,k,t])], shape1 = a, shape2 = b,
                                                                                 MMa = MMa, MMb = MMb, sigma=sigma[t], log = TRUE)))
    },
    error = function(e){sum = Inf})
    return(sum)
  }
  return(f)  
}

.NLL_ExpMMtailL_scale = function(x, z_nkt, pi, times, MMa, MMb, sigma, shape) 
{
  f = function(scale) 
  { 
    n_t = length(times)
    sum = 0
    if ((scale<=0)|(scale>0.8)){ return(Inf) } 
    
    for (t in 1:n_t){
      tryCatch({
        sum = sum - (z_nkt[,1,t][!is.na(z_nkt[,1,t])] %*% ( log(pi[[1]][t]) +
                                                              d_expMMtail_convolved(x[[t]][!is.na(z_nkt[,1,t])], shape = shape[t], 
                                                                                    scale = scale, MMa = MMa, MMb = MMb, sigma=sigma[t], log = TRUE)))
      },
      error = function(e){sum = Inf})
    }
    return(sum)
  }
  return(f)  
}

.NLL_ExpMMtailL_shape = function(x, z_nkt, pi, MMa, MMb, sigma, scale, t) # pass x[[t]] as x here
{
  f = function(shape) 
  { 
    if ((shape<=1e-2)|(shape>15)){ return(Inf) } 
    
    tryCatch({
      sum = - (z_nkt[,1,t][!is.na(z_nkt[,1,t])] %*% ( log(pi[[1]][t]) +
                                                        d_expMMtail_convolved(x[!is.na(z_nkt[,1,t])], shape = shape, 
                                                                              scale = scale, MMa = MMa, MMb = MMb, sigma=sigma[t], log = TRUE)))
    },
    error = function(e){sum = Inf})
    return(sum)
  }
  return(f)  
}

.NLL_ExpMMtailG_shape = function(x, z_nkt, pi, MMa, MMb, sigma, scale, t) # pass x[[t]] as x here
{ # for inferring shape parameters using fixed scale
  f = function(y) 
  { 
    if ((y[1]<=1e-2)|(y[1]>15)){ return(Inf) } 
    
    tryCatch({
      sum = - (z_nkt[,1,t][!is.na(z_nkt[,1,t])] %*% ( log(pi[[1]][t]) +
                                                        d_expMMtail_convolved(x[!is.na(z_nkt[,1,t])], shape = y[1], 
                                                                              scale = scale, MMa = MMa, MMb = MMb, sigma=sigma[t], log = TRUE)))
    },
    error = function(e){sum = Inf})
    return(sum)
  }
  
  return(f)  
}

.NLL_ExpMMtailG_scale = function(x, z_nkt, pi, times, MMa, MMb, sigma, shape) 
{ # for inferring scale parameter that is shared between timepoints with fixed shape
  f = function(y) # params: scale
  { 
    n_t = length(times)
    sum = 0
    
    if ((y[1]<=0)|(y[1]>0.8)){ return(Inf) } 
    
    for (t in 1:n_t){
      tryCatch({
        sum = sum - (z_nkt[,1,t][!is.na(z_nkt[,1,t])] %*% ( log(pi[[1]][t]) +
                                                              d_expMMtail_convolved(x[[t]][!is.na(z_nkt[,1,t])], shape = shape[t], 
                                                                                    scale = y[1], MMa = MMa, MMb = MMb, sigma=sigma[t], log = TRUE)))
      },
      error = function(e){sum = Inf})
    }
    return(sum)
  }
  
  return(f)  
}

.NLL_sigma = function(x, fit, tail, t, transformed) # pass x[[t]] as x here
{
  f = function(y) # NLL
  {
    shape = unlist(fit$shape)  
    scale = unlist(fit$scale)
    
    if (any(y <= 0) | any(y > 0.5)){ return(Inf) } 
    
    sum = 0
    tryCatch({
      for (k in 1:length(fit$pi)){
        
        if (k!=1){
          ak = fit$a[[k-1]] 
          bk = fit$b[[k-1]]  
        }
        
        if (transformed){
          if (k==1){
            if (!tail){density=0} else {density = d_expMMtail_convolved(x[!is.na(fit$z_nkt[,k,t])], shape = shape[t], scale = scale, 
                                                                        MMa = fit$MM_scale$MM_scale_a, MMb = fit$MM_scale$MM_scale_b, 
                                                                        sigma = y, log = TRUE)}}
          else
          {density = d_expMMbeta_convolved(x[!is.na(fit$z_nkt[,k,t])], shape1 = ak[t], shape2 = bk[t],
                                           MMa = fit$MM_scale$MM_scale_a, MMb = fit$MM_scale$MM_scale_b,
                                           sigma=y, log = TRUE)} 
        } else {
          if (k==1){
            if (!tail){density=0} else {density = d_pareto_convolved(x[!is.na(fit$z_nkt[,k,t])],  
                                                                     shape = shape[t], scale = scale, sigma = y, log = TRUE)}
          } 
          else
          {density = d_beta_convolved(x[!is.na(fit$z_nkt[,k,t])], a = ak[t], b = bk[t],
                                      sigma = y, log = TRUE)}
          
        }
        if ((k==1) && (!tail)){term = 0} else {term = fit$z_nkt[,k,t][!is.na(fit$z_nkt[,k,t])] %*% (log(fit$pi[[k]][t]) + density)}
        sum = sum - term
      }
    },
    error = function(e){sum = Inf})
    return(sum)
  } 
  return(f)
}