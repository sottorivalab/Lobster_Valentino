# Return latent variables (re-computed), plus other variables.
latent_vars = function(x, clusters_tibble = NULL, log_MMscaled = FALSE, MM_scale_a = NULL, MM_scale_b = NULL, untf = FALSE) {

  if(all(!is.null(clusters_tibble)))
    x$Clusters = clusters_tibble

  # Extract Beta values
  B = mobster:::.params_Beta(x)
  a = B$a
  b = B$b
  
  if (nrow(B)>0){
    names(a) = names(b) = B$cluster
  }

  # Extract Tail values
  shape = ifelse(x$fit.tail, mobster:::.params_Pareto(x)$Shape, NA)
  scale = ifelse(x$fit.tail, mobster:::.params_Pareto(x)$Scale, NA)

  sigma = x$sigma
  
  # Extract Mixing Proportions
  pi = mobster:::.params_Pi(x)

  # Reconstruct latent variables, responsibilities and weighted densities
  N = nrow(x$data)
  z_nkt      = array(0, c(N,x$K,x$n_t))  
  #colnames(z_nkt) = names(pi)
  
  pdf.w = z_nkt

  # Get log density per component
  for (k in 1:x$K){
    for (t in 1:x$n_t){
      pdf.w[,k,t] = mobster::ddbpmm(x,
                                        data = x$data$VAF,
                                        components = paste(k),
                                        a = a,
                                        b = b,
                                        pi = pi,
                                        shape = shape,
                                        scale = scale,
                                        log = TRUE,
                                        log_MMscaled = log_MMscaled,
                                        MM_scale_a = MM_scale_a,
                                        MM_scale_b = MM_scale_b,
                                        sigma = sigma[t],
                                        untf = untf,
                                        time_point = t 
      )
    }
  }
  
  Z       = apply(pdf.w, c(1,3), mobster:::.log_sum_exp)
  z_nkt   = sweep(pdf.w,MARGIN=c(1,3),Z,FUN="-")  # fit$pdf.w - Z 
  z_nkt   = apply(z_nkt, c(2,3), exp) # Exponentiate to get actual probabilities
  NLL     = -sum(Z, na.rm=TRUE)  # Evaluate the NLL (ignores NA)
  
  NLL.hard = calc_NLL.hard(pdf.w, x) 
  
  return(list(NLL = NLL, NLL.hard = NLL.hard, z_nkt = z_nkt, pdf.w = pdf.w, a = a, b = b, pi = pi, shape = shape, scale = scale, sigma = sigma))
}

# Compute hard clustering assignments
latent_vars_hard_assignments = function(lv) {
  # Cluster labels of each data point. Each data point is assigned to the cluster
  # with the highest posterior responsibility. Hard assignment.
  
  resp_prods = apply(lv$z_nkt, c(1,2), function(x){mean(x,na.rm = TRUE)}) # was prod # multiply responsibilities across timepoints (CHECK!)
  unlist(apply(resp_prods, 1, function(x) {names(lv$pi)[which(x == max(x, na.rm = TRUE))[1]]})) 
}

z_nkt_hard_clustered_mask = function(z_nkt, fit){
  # create mask from hard-clustering selection to add to z_nkt
  resp_prods = apply(z_nkt, c(1,2), function(x){mean(x,na.rm = TRUE)}) # was prod # multiply responsibilities across timepoints (CHECK!)
  chosen_k = apply(resp_prods, 1, function(x) {which(x == max(x, na.rm = TRUE))[1]})
  z_nkt_mask_list = lapply(chosen_k, function(x){my_arr = array(-Inf, dim = c(1, fit$K, fit$n_t)); my_arr[,x,] = 0; my_arr})
  z_nkt_mask = do.call(abind, c(z_nkt_mask_list, along = 1))
  
  return(z_nkt_mask)
}

calc_NLL.hard = function(pdf.w, fit){
  # calculate 
  mask       = z_nkt_hard_clustered_mask(pdf.w, fit)
  Z.hard     = apply((pdf.w + mask), c(1,3), mobster:::.log_sum_exp)
  NLL.hard   = -sum(Z.hard, na.rm=TRUE)  # Evaluate the NLL (ignores NA)
  
  return(NLL.hard)
}

# Compute scores for model selection 
latent_vars_scores = function(lv, K, n_t, tail, cluster, convolved)
{
  ##==============================
  # Scores for model selection   #
  ##==============================
  
  # K = Kbeta + 1
  pi_components = ((K-1) + as.numeric(tail))*n_t
  Beta_components = ((K-1) * 2)*n_t
  Pareto_components = (as.numeric(tail))*(n_t + 1)
  Sigma_components = ifelse(convolved,n_t,0)
  
  numParams = pi_components + Beta_components + Pareto_components + Sigma_components
  
  N = nrow(lv$z_nkt)
  
  BIC <-
    2 * lv$NLL + numParams * log(N)          # BIC = -2*ln(L) + params*ln(N)
  BIC.hard1 <-
    2 * lv$NLL.hard + numParams * log(N)     # BIC with hard-clustering mask imposed
  BIC.hard2 <-
    2 * lv$NLL.hard + (numParams * log(N)/K)     # BIC with hard-clustering mask imposed
  AIC <-
    2 * lv$NLL + 2 * numParams               # AIC = -2*ln(L) + 2*params
  AIC.hard <-
    2 * lv$NLL.hard + 2 * numParams               # AIC = -2*ln(L) + 2*params
  
  # Integrated Complete Likelihood criterion -- uses standard entropy
  entropy <- -sum(lv$z_nkt * log(lv$z_nkt), na.rm = TRUE)
  ICL <- BIC + entropy
  
  # Integrated Complete Likelihood criterion with reduced entropy (only for
  # latent variable that involve subclones -- i.e., Betas). I think this is also a sort of
  # conditional entropy where we condition on the MAP estimate of a mutation being part
  # of a clone, rather than tail.
  
  # Take the MAP estimates of z_nk, and select only entries that are assigned
  # to a Beta component (i.e. those with arg_max != Tail)
  if (K>1){
    cz_nkt = lv$z_nkt[cluster != 'Tail', 2:ncol(lv$z_nkt),,drop = FALSE]
    
    # This is un-normalized -- we compute the empirical normalizing constant (C)
    C   = apply(cz_nkt, c(1,3), sum) # if any responsibilities across components are NA, isn't counted
    cz_nkt <- sweep(cz_nkt, 2, C, "/", check.margin=F) #cz_nk = cz_nk/C
    
    # The reduced entropy is the entropy of this distribution
    rentropy = -sum(cz_nkt  * log(cz_nkt), na.rm = TRUE)
    
    # Integrated Complete Likelihood criterion with reduced entropy
    reICL <- BIC + rentropy
  } else {
    cz_nkt = C = rentropy = reICL = NA
  }
  
  scores = data.frame(
    NLL = lv$NLL,
    NLL.hard = lv$NLL.hard,
    BIC = BIC,
    BIC.hard1 = BIC.hard1,
    BIC.hard2 = BIC.hard2,
    AIC = AIC,
    AIC.hard = AIC.hard,
    entropy = entropy,
    ICL = ICL,
    reduced.entropy = rentropy,
    reICL = reICL,
    size = numParams
  )
  
  scores
}

