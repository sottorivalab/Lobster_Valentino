# Return initialization parameters for the fit
.initializer = function(X, K, n_t, tail, init, pareto.shape = list(min.val = 0.01, max.val = 5),
                        MM_scale_a = NULL, MM_scale_b = NULL, log_MMscaled = NULL, 
                        untf_xstar = NULL, xstar = NULL, times = c(0))
{
  
  if(is.list(init)) {

    stopifnot(c('beta', 'shape', 'scale') %in% names(init))
    stopifnot(c('a', 'b') %in% rownames(init$beta))

    return(init)
  }

  # Sampler for Beta parameters with random mean.
  # Reject samples if a/ b are negative
  bsamp = function(m = runif(1)) {
    repeat{
      v = runif(n = 1, min = 0.0001, max = 0.01)
      p = c(mobster:::.estBetaParams(m, v), mean=m, var=v)
      if(all(p > 0)) return(p)
    }
  }

  init = 'random' # setting for devel purposes
  if(is.character(init) && init == 'random'){
    beta = sapply(runif(K*n_t), bsamp) %>% unlist() %>% array(c(4,K,n_t))
  }  

  # Prepare tibbles for Betas
  Betas = NULL

  for(i in 1:K)
    Betas =  dplyr::bind_rows(
      Betas,
      tibble::tribble(
        ~cluster, ~type, ~fit.value, ~init.value,
        paste0('C',i), "a", rep(NA,n_t), unlist(beta[1,i,]) %>% as.vector,
        paste0('C',i), "b", rep(NA,n_t), unlist(beta[2,i,]) %>% as.vector,
        paste0('C',i), "Mean", rep(NA,n_t), unlist(beta[3,i,]) %>% as.vector,
        paste0('C',i), "Variance", rep(NA,n_t), unlist(beta[4,i,]) %>% as.vector
      )
    )

  # Tail, sample and make a tibble
  Pareto = NULL

  if(tail)
  {
    shape = runif(n_t, min = pareto.shape$min.val, max = pareto.shape$max.val)
    if (!is.null(xstar)){
      scale = xstar
    } else {
      scale = min(X,na.rm = TRUE) - 1e-9
    }
    mv = mobster:::.MeanVarPareto(shape, scale)
    Pareto = tibble::tribble(
      ~cluster, ~type, ~fit.value, ~init.value,
      'Tail', "Shape", rep(NA,n_t), shape,
      'Tail', "Scale", NA, unlist(scale) %>% as.vector,
      'Tail', "Mean", rep(NA,n_t), unlist(mv$mean) %>% as.vector,
      'Tail', "Variance", rep(NA,n_t), unlist(mv$var) %>% as.vector
    )
  }

  # Mixing proportions -- uniform distribution
  K = K + as.numeric(tail)

  Clusters = dplyr::bind_rows(Pareto, Betas)
  
  # pi will not be treated as coefficients, but rather as simple values, for practicality
  # trend can be calculated later through regression
  pi.Pareto = tibble::tribble(
    ~cluster, ~type, ~fit.value, ~init.value,
    'Tail', "Mixing proportion", rep(NA,n_t), ifelse(rep(tail,n_t), rep(1/K,n_t), rep(0,n_t))
  )
  
  pi.Betas = NULL
  for(x in unique(Betas$cluster))
    pi.Betas = dplyr::bind_rows(pi.Betas,
                                tibble::tribble(
                                  ~cluster, ~type, ~fit.value, ~init.value,
                                  x, "Mixing proportion", rep(NA,n_t), rep(1/K,n_t)
                                ))

  dplyr::bind_rows(Pareto, Betas, pi.Pareto, pi.Betas)
}

# initial hard clustering through kmeans
.kmeans_resp_init = function(X, Kbeta, n_t, tail, z_nkt){
  
  kmeans_NA_handler = function(data, ...){
    # handles NA in data while preserving vectorisation
    NA_mask = is.na(data)
    nonNA_idx = which(!NA_mask)
    term = array(NA, c(length(data)))
    filtered_data = data[!NA_mask]
    
    filtered_kmeans = kmeans(filtered_data,...) 
    
    return(list(nonNA_idx = nonNA_idx, 
                cluster = filtered_kmeans$cluster, 
                centers = filtered_kmeans$centers))
  }

  # assign initial hard clusters through k-means clustering
  # cluster with lowest mean assumed to be tail if fit.tail==TRUE
  Kcl = Kbeta + ifelse(tail,1,0)
  znkt_cl = array(data = NA, dim = dim(z_nkt))
  if (Kcl > 1){
    for (t in 1:n_t){
      # Cluster their x-coordinates
      clus = kmeans_NA_handler(X[[t]], Kcl, nstart=10)
      cnt = as.vector(clus$centers)
      cl = as.vector(clus$cluster)
      kord = order(cnt) # order returned clusters by ascending mean
     # new_cl = kord[cl] # clusters in our order
      new_cl = match(cl,kord)
      if (tail==FALSE){new_cl = new_cl + 1} # offset to skip tail
      
      # set to 1 at selected cluster
      for (i in 1:length(as.vector(clus$nonNA_idx))){
        znkt_cl[as.vector(clus$nonNA_idx)[i],new_cl[i],t] = 1
      }
      non_empty_rows <- apply(!is.na(znkt_cl[,,t]), 1, any)
      # Fill NAs with zero in non-empty rows
      znkt_cl[non_empty_rows,,t] <- apply(znkt_cl[non_empty_rows,,t], 2, function(x) replace(x,is.na(x),0))
    }
  } else {znkt_cl = z_nkt}
  
  return(znkt_cl) 
}


# Decide if the fit stops or not
.stoppingCriterion = function(i, prevNLL, NLL, prevpi, pi, fit.type, epsilon, isDebug, K)
{
  # Compute NLL difference after ith iteration
  NLL.Diff  <- prevNLL - NLL

  # Fit by MLE
  if(fit.type == 'MLE')
  {
    # Correctness: we know by EM properties that the LL should always increase
    if (NLL.Diff < -1e-10)
    {
      stop("Negative log likelihood increases: ", NLL.Diff, ".\nSomething is wrong, aborting!")
    }

    # Check for convergence
    stopping = NLL.Diff < epsilon
    
  }

  # Fit by MM
  if(fit.type == 'MM')
  {
    # With MM in general we are not interested in the logLik
    # and we only look for variations in the actual mixing
    # proportions. When the fit is however single-cluster, 
    # we need to check also the logLik because the proportions
    # never change by definition (do the abs beacuse Jensen's ineq. is invalid)
    
    
    if(K > 1)
    {
      # Compute pi's difference after ith iteration
      #pi.Diff  <- abs(prevpi - pi)
      pi.Diff <-  mapply('-', prevpi, pi, SIMPLIFY = FALSE) %>% lapply(abs) %>% unlist()
  
      # Check for convergence -- all pi's have changed less than epsilon
      stopping = all(pi.Diff < epsilon)
    }
    else
      stopping = abs(NLL.Diff) < epsilon
  }

  # Some printing
  # NEEDS CHANGING
 # if (isDebug)
 # {
 #   flush.console()
 #   pi = round(pi, 3)
 #   NLL = round(NLL, 2)
 #   NLL.Diff = format(NLL.Diff, scientific = T)
 #   pi.Diff = max(prevpi - pi)
#
 # }

  return(stopping)
}

# Peaks detection function
.find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

# Use the log sum exp trick for having numeric stability
.log_sum_exp <- function(x) {
  # Computes log(sum(exp(x))
  offset <- max(x)
  s <- log(sum(exp(x - offset))) + offset
  i <- which(!is.finite(s))
  if (length(i) > 0) { s[i] <- offset }
  return(s)
}

# Switch from mean/ var to a and b for Beta
.estBetaParams <- function(mu, var)
{
  if(var == 0) var = 1e-9
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(a = alpha, b = beta))
}

# Switch from a and b to mean/ var for Beta
.MeanVarBeta = function(a, b) {
  m = a/(a+b)
  return(list(mean = m, var = (m*(1-m))/(a+b+1)))
}

# Switch from shape and scale to mean/ var for Pareto
.MeanVarPareto = function(shape, scale) {
  # need to calc true values (from poly coeffs) before passing to this function
  meanPareto = ifelse(shape > 1, (shape * scale) / (shape - 1), Inf)
  varPareto = ifelse(shape > 2, (shape * scale**2) / ((shape - 1)**2 * (shape - 2)), Inf)
  
  return(list(mean = meanPareto, var = varPareto))
}

# Extract Beta parameters
.params_Beta = function(x, init = FALSE)
{
  if(init) x$Clusters$fit.value = x$Clusters$init.value
  
  x$Clusters %>%
    dplyr::filter(cluster != 'Tail', type == 'a' | type == 'b') %>%
    dplyr::select(-init.value) %>%
    tidyr::spread(key = type, value = fit.value)
}

# Extract Pareto parameters
.params_Pareto = function(x, init = FALSE)
{
  if(init) x$Clusters$fit.value = x$Clusters$init.value
  
  suppressWarnings(
    x$Clusters %>%
      dplyr::filter(cluster == 'Tail', type == 'Shape' | type == 'Scale') %>%
      dplyr::select(-init.value) %>%
      tidyr::spread(key = type, value = fit.value)
  )
}

# Extract mixing proportions parameters
.params_Pi = function(x, init = FALSE)
{
  if(init) x$Clusters$fit.value = x$Clusters$init.value
  
  v = x$Clusters %>%
    dplyr::filter(type == 'Mixing proportion') %>%
    dplyr::select(-init.value) %>%
    tidyr::spread(key = type, value = fit.value)
  
  # pi = v$`Mixing proportion`
  # names(pi) = v$cluster
  
  pi = pio:::nmfy(v$cluster, v$`Mixing proportion`)
  
  # ord.pi = c(pi['Tail'], pi[names(pi) != 'Tail'])
  ord.pi = c(
    pi['Tail'],
    pi[sort(names(pi)[names(pi) != 'Tail'], decreasing = FALSE)]
    )
  
  ord.pi
}

# Set Beta parameters
.set_params_Beta = function(fit, a, b, times)
{
  names.BetaC = names(a)
  if(any(is.null(names.BetaC))) stop("params -- named vector required?")
  
  fit$Clusters = fit$Clusters %>%
    dplyr::mutate(
      fit.value = 
        ifelse(
          (cluster %in% names.BetaC) & (type == 'a'), 
          a[cluster],
          fit.value
        )
    )
  
  fit$Clusters = fit$Clusters %>%
    dplyr::mutate(
      fit.value = 
        ifelse(
          (cluster %in% names.BetaC) & (type == 'b'), 
          b[cluster],
          fit.value
        )
    )
  
  for (s in names.BetaC)
  {
    
    for (i in 1:length(fit$times)){
      if (i==1){
        mv = as.data.frame(mobster:::.MeanVarBeta(a[[s]][i], b[[s]][i]))
        }
      else {
        new = as.data.frame(mobster:::.MeanVarBeta(a[[s]][i], b[[s]][i]))
        mv = bind_rows(mv, new)
        }
      
    }
    
    fit$Clusters = fit$Clusters %>%
      dplyr::mutate(
        fit.value = 
          ifelse(
            (cluster == s) & (type == 'Mean'), 
            list(mv$mean),
            fit.value
          )
      )
    
    fit$Clusters = fit$Clusters %>%
      dplyr::mutate(
        fit.value = 
          ifelse(
            (cluster %in% s) & (type == 'Variance'), 
            list(mv$var),
            fit.value
          )
      )
  }
  
  fit
}

# Set Pareto parameters
.set_params_Pareto = function(fit, shape, scale, times)
{
  if(!fit$fit.tail) return(fit)
  
 #  commented code causes program to crash and I don't know why :|
 # mv = mobster:::.MeanVarPareto(shape_true, scale)
  
  fit$Clusters[
    fit$Clusters$cluster == 'Tail' & fit$Clusters$type == 'Shape',
    'fit.value'
    ] = tibble(shape)
  
  fit$Clusters[
    fit$Clusters$cluster == 'Tail' & fit$Clusters$type == 'Scale',
    'fit.value'
    ] =  list(scale)
  
  fit
}

# Set Exponentiated MM-scaled Beta parameters
.set_params_ExpMMbeta = function(fit, a, b)
{
  names.BetaC = names(a)
  if(any(is.null(names.BetaC))) stop("params -- named vector required?")
  
  fit$Clusters = fit$Clusters %>%
    dplyr::mutate(
      fit.value = 
        ifelse(
          (cluster %in% names.BetaC) & (type == 'a'), 
          a[cluster],
          fit.value
        )
    )
  
  fit$Clusters = fit$Clusters %>%
    dplyr::mutate(
      fit.value = 
        ifelse(
          (cluster %in% names.BetaC) & (type == 'b'), 
          b[cluster],
          fit.value
        )
    )
  
  return(fit)
}

# Set  parameters
.set_params_ExpMMtail = function(fit, shape, scale)
{
  if(!fit$fit.tail) return(fit)
  
  fit$Clusters[
    fit$Clusters$cluster == 'Tail' & fit$Clusters$type == 'Shape',
    'fit.value'
    ] = tibble(as.list(shape))
  
  fit$Clusters[
    fit$Clusters$cluster == 'Tail' & fit$Clusters$type == 'Scale',
    'fit.value'
    ] =  list(as.list(scale))
  
  fit
}

# Set mixing proportions parameters
.set_params_Pi = function(fit, pi)
{
  fit$Clusters = fit$Clusters %>%
    dplyr::mutate(
      fit.value = ifelse(
        type == 'Mixing proportion', 
        pi[cluster],
        fit.value
      )
    )
  
  fit
}

# Compute SSE of a fit
.compute_fit_sqerr = function(x, 
                              binning = 1e-2, 
                              log_MMscaled = FALSE, 
                              MM_scale_a = NULL, 
                              MM_scale_b = NULL,
                              untf = FALSE)
{
  densities = mobster:::template_density(
    x,
    x.axis = seq(binning, 1-binning, by = binning), # Restricted for numerical errors
    binwidth = binning,
    reduce = TRUE,
    log_MMscaled = log_MMscaled,
    MM_scale_a = MM_scale_a,
    MM_scale_b = MM_scale_b,
    untf = untf)
  
  densities = tibble::as_tibble(densities)
  densities = densities %>% group_by(x) %>% summarise(y = sum(y), cluster = 'f(x)')
  
  # Empirical density
  empirical = hist(x$data$VAF, breaks = seq(0, 1, binning), plot = FALSE)$density
  empirical = empirical[-length(empirical)]
  empirical = empirical * binning # adjust for binwidth
  
  # Error
  error = densities
  error$cluster = 'e(x)'
  error$y = (densities$y - empirical)^2
  error$cum.y = cumsum(error$y)
  
  error
}

.get_clusters_labels = function(x)
{
  unique(x$data$cluster)
}


# Very fast analysis, raw results but good to get a grasp of the data
# - max 2 clones (1 subclone);
# - random start
# - 2 samples
# - mild epsilon/ maxIter
# - not parallel
template_parameters_fast_setup = function()
{
  return(
    list(
      K = 1:2,
      samples = 2,
      init = 'random',
      tail = c(TRUE, FALSE),
      epsilon = 1e-6,
      maxIter = 100,
      fit.type = 'MM',
      seed = 12345,
      model.selection = 'reICL',
      trace = FALSE,
      parallel = FALSE,
      pi_cutoff = 0.02,
      N_cutoff = 10
    )
  )
}

auto_setup = function(x)
{
  if(x == "FAST")  return(template_parameters_fast_setup())
  
  stop("Auto setup unknown: use \"FAST\".")
}

.update_fit_table = function(fit, log_MMscaled=FALSE, MM_scale_b=NULL)
{
    # Update fit table
    fit = mobster:::.set_params_Beta(fit, fit$a, fit$b, fit$times)
    fit = mobster:::.set_params_Pareto(fit, fit$shape, fit$scale, fit$times)
    fit = mobster:::.set_params_Pi(fit, fit$pi)
    
    return(fit)
}    