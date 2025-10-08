# The actual fitting function
.dbpmm.EM <-
  function(X,
           K = 3,
           init = 'peaks',
           tail = TRUE,
           epsilon = 1e-10,
           maxIter = 1000,
           is_verbose = FALSE,
           fit.type = 'MM',
           trace = FALSE,
           pi_cutoff = 0.02,
           N_cutoff = 10,
           log_MMscaled = FALSE,
           MM_scale_a = NULL,
           MM_scale_b = NULL,
           xstar = NULL,
           seed = 12345,
           genSA_steps = 2,
           description = "My MOBSTER model.",
           Gopt = 'DE',
           untf_convolved = FALSE)
  {
     set.seed(seed)
     skip_calc = FALSE
     
      if(log_MMscaled && !(tail==FALSE) && is.null(xstar)){
        warning('xstar not specified and is therefore taken to be min(VAF) - this may lead to problems fitting transformed tail\n')
      }
  
    stopifnot(fit.type %in% c('MLE', 'MM'))
    stopifnot(tail | K > 0)
    
    .onLoad(NULL, NULL)
    
    # suppressMessages(require(tidyverse))
    # suppressMessages(require(pio))
    # suppressMessages(require(crayon))

    ##=============================##
    # Create a BetaParetoMM object  #
    ##=============================##
    fit           = list()
    class(fit) <- "dbpmm"

    fit$data      = X
    fit$Call      = match.call()
    fit$description = description

    fit$fit.type  = fit.type # MLE or MM

    fit$fit.tail  = tail
    fit$Kbeta     = K
    fit$K         = K + 1 # One extra pareto Mode -- specially assigned to slot #1

    fit$N         = nrow(X)   # Number of samples
    fit$n_t       = get_no_timepoints(X)
    fit$times     = fit$data$times[1,]
    fit$z_nkt     = array(0, c(fit$N,fit$K,fit$n_t))        # Responsibilities
    fit$N.k       = NULL       # Clustering assignments

    fit$pdf.w     = array(0, c(fit$N,fit$K,fit$n_t))        # Weighted PDFs
    fit$all.NLL   = vector(mode = "numeric")           # Hold NLL for all EM iterations
    fit$NLL       = 1e+40                            # Initialize Negative Log Likelihood

    fit$labels    = NULL # hard clustering assignments
    fit$trace     = NULL # trace for the fit

    # transform xstar to log_MM transformed space
    if ((log_MMscaled)&&(!is.null(xstar))){
      xstar_t = MM_scale_a*log(xstar) + MM_scale_b
    } else {xstar_t = xstar}
    
    # Names of components
    names.BetaC = paste('C', 1:K, sep = '')
    names.ParetoC = 'Tail'

    if (log_MMscaled){
      fit$MM_scale = list(MM_scale_a = MM_scale_a, MM_scale_b = MM_scale_b)
    }
    
    # Compute initial conditions
    fit$Clusters  = mobster:::.initializer(X$VAF, 
                                           K = fit$Kbeta, 
                                           n_t = fit$n_t,
                                           tail = tail, 
                                           init = init, 
                                           MM_scale_a = MM_scale_a, 
                                           MM_scale_b = MM_scale_b, 
                                           log_MMscaled = log_MMscaled,
                                           untf_xstar = xstar,
                                           xstar = xstar_t,
                                           times = fit$times)
    fit$Clusters$fit.value = fit$Clusters$init.value
    
    
    # Initialise Gaussian std as equal to untransformed xstar
    fit$sigma = rep(xstar,fit$n_t)

    # Extract Beta values
    fit$a = mobster:::.params_Beta(fit)$a
    fit$b = mobster:::.params_Beta(fit)$b

    names(fit$a) = names(fit$b) = names.BetaC

    # Extract Tail values
    fit$shape = rep(NA,fit$n_t)
    fit$scale = NA
    if(fit$fit.tail)
    {
      fit$shape = mobster:::.params_Pareto(fit)$Shape
      fit$scale = mobster:::.params_Pareto(fit)$Scale
    }

    # Extract Mixing Prop
    fit$pi = mobster:::.params_Pi(fit)
    
    # saving kmeans hard-clustering for plotting
    fit$z_nkt_kmeans = fit$z_nkt
    fit$data$kmeans_cluster = mobster:::.cluster_kmeans(fit)
    
    fit$scores = NULL
    
    # create two fit objects, one initialised with log-scaled kmeans and one with normal kmeans
    # after one EM iteration, compare NLL and choose better to continue with
    fits = list()
    fits[[1]] = fits[[2]] = fit # initialise with same fit
    
    # log-scaled kmeans
    # transform VAF in this function only
    fits[[1]]$z_nkt = mobster:::.kmeans_resp_init(mobster:::.log_transform(X)$VAF_ln_scaled, 
                                          K = fit$Kbeta, 
                                          n_t = fit$n_t,
                                          tail = tail,
                                          z_nkt = fit$z_nkt)
    
    # unscaled kmeans
    fits[[2]]$z_nkt = mobster:::.kmeans_resp_init(X$VAF, 
                                          K = fit$Kbeta, 
                                          n_t = fit$n_t,
                                          tail = tail,
                                          z_nkt = fit$z_nkt)
    
    ##=========================================
    # Run Expectation Maximization  algorithm #
    ##=========================================
    stopping = FALSE
    
    i = 0
    while ((!stopping)&(i<=maxIter)){
        print(i)
      
        if (i==0){
          fit = fits[[1]]  
          fit$kmeans_type = 'scaled'
        } else if (i==1){
          fit = fits[[2]]
          fit$kmeans_type = 'unscaled'
        }
          
        # Convergence:
        # - NLL with MLE
        prevNLL = fit$NLL
        prevpi = fit$pi
        fit = EM_step(fit, X, i, 
                    log_MMscaled,
                    MM_scale_a,
                    MM_scale_b,
                    untf_convolved,
                    names.BetaC,
                    names.ParetoC,
                    genSA_steps,
                    Gopt, 
                    trace, 
                    tail)
        
        # an E-step for evaluation of NLL
        for (k in 1:fit$K){
          for (t in 1:fit$n_t){
            fit$pdf.w[,k,t] = mobster::ddbpmm(fit,
                                    data = fit$data$VAF,
                                    components = k,
                                    a = fit$a,
                                    b = fit$b,
                                    pi = fit$pi,
                                    shape = fit$shape,
                                    scale = fit$scale,
                                    log = TRUE,
                                    log_MMscaled = log_MMscaled,
                                    MM_scale_a = MM_scale_a,
                                    MM_scale_b = MM_scale_b,
                                    sigma = fit$sigma[t],
                                    time_point = t,
                                    untf = untf_convolved
                                    )
          }
        }
        
        # for calculating pi (for each t)
        Z           = apply(fit$pdf.w, c(1,3), mobster:::.log_sum_exp)
        fit$z_nkt   = sweep(fit$pdf.w,MARGIN=c(1,3),Z,FUN="-")  # fit$pdf.w - Z 
        fit$z_nkt   = apply(fit$z_nkt, c(2,3), exp) # Exponentiate to get actual probabilities
        if (any(is.infinite(fit$z_nkt)))
          stop('Error? Some latent variables (z_nkt) are Infinite.')
        fit$NLL    = -sum(Z, na.rm=TRUE)  # Evaluate the NLL (ignores NA)
      
      
        if (i==0){
          fits[[1]] = fit
        } else if (i==1){
          fits[[2]] = fit
        }
        
        ## Convergency test 
                    
        if (i==1){
          # choose best fit out of fits[[1]] and fits[[2]]
          if (fits[[2]]$NLL > fits[[1]]$NLL){
            fit = fits[[1]]
          } else {
            fit = fits[[2]]
          }
          
        } else if (i>1) {
          stopping=mobster:::.stoppingCriterion(i,
                                 prevNLL,
                                 fit$NLL,
                                 prevpi,
                                 fit$pi,
                                 fit.type,
                                 epsilon,
                                 is_verbose,
                                 fit$K)
        }
        
        i = i + 1
    }
    #End of Expectation Maximization loop.
    if(is_verbose) print("EM completed!")
    
    ########### Add names to the estimated variables for clarity, populate summary matrices
    names(fit$pi) =  colnames(fit$z_nkt) = colnames(fit$pdf.w) = c(names.ParetoC, names.BetaC)
    names(fit$a) = names(fit$b) = names.BetaC
   
    fit = mobster:::.update_fit_table(fit,log_MMscaled,MM_scale_b) # update Clusters object
     
    # finds best permutation of beta clusters per timepoint, permutes z_nkt, pdf.w and parameters of fit
    # uses z_nkt from last E step - should prob change to do another E step to be completely up to date
    if (fit$Kbeta > 1){
      fit = mobster:::rename_Beta_clusters_permutation(fit)
    }
    fit = mobster:::.update_fit_table(fit,log_MMscaled,MM_scale_b) 

    # Cluster labels of each data point. Each data point is assigned to the cluster
    # with the highest posterior responsibility. Hard assignment.
    resp_prods = apply(fit$z_nkt, c(1,2), function(x){mean(x,na.rm = TRUE)})#prod(x,na.rm = TRUE)}) # multiply responsibilities across timepoints (CHECK!)
    fit$data$cluster =  unlist(apply(resp_prods, 1,
                                     function(x) {names(fit$pi)[which(x == max(x, na.rm = TRUE))[1]]})) 
    
    # Summary numbers
    fit$N.k = rep(0, fit$K)
    names(fit$N.k) = names(fit$pi)
    obFreq = table(fit$data$cluster)

    fit$N.k[names(obFreq)] = obFreq

    ##==============================
    # Scores for model selection   #
    ##==============================
    
    fit$scores = latent_vars_scores(
      latent_vars(fit, # permutation should not affect scores
                  log_MMscaled = log_MMscaled,
                  MM_scale_a = MM_scale_a,
                  MM_scale_b = MM_scale_b), # Extract latent variables
      fit$K,
      fit$n_t,
      fit$fit.tail,
      fit$data$cluster,
      (log_MMscaled|untf_convolved) # whether convolved distributions
      )
    
    if (is_verbose) {
      print.dbpmm(fit)
      cat('\n')
    }
    
    # Now apply the cluster-selection heuristic in function choose_clusters
    fit = get_cluster_mean_and_var(fit, log_MMscaled) 
    
    fit = choose_clusters(fit,  pi_cutoff = pi_cutoff, N_cutoff = N_cutoff,
                  log_MMscaled = log_MMscaled,
                  untf_convolved = untf_convolved,
                  MM_scale_a = MM_scale_a,
                  MM_scale_b = MM_scale_b)

    fit = get_cluster_mean_and_var(fit, log_MMscaled)
    # ... and re-order the Beta cluster ID by mean ...
    #fit = rename_Beta_clusters(fit)
    fit = hard_clustered_VAF_CCF(fit)
    
    return(fit)
  }


