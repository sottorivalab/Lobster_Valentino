
EM_step <- function(fit, 
                    X,
                    i, 
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
  {
  
    ##################
#
    if ((i>1)|(fit$Kbeta==1)){
      
        ##################
      # Store trace to visualize fit
      # current =  data.frame(step = i, NLL = fit$NLL,
 
      if(trace){
        step.density = fit$Clusters
        step.density$step = i
        fit$trace = dplyr::bind_rows(fit$trace, step.density)
      }
 
      ##===================
      #       E-Step      #
      ##===================
      # When pi(Pareto) --> 0 the MLE fit for shape --> 0/0 = NaN and thus at some point we get NaN here
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
    } else { fit$NLL = 1e+20 } 
    
    fit$all.NLL   <-
      c(fit$all.NLL, fit$NLL)    # Keep all NLL in a vector
    
    ##===================
    #       M-Step      #
    ##===================
    N.kt   <- apply(fit$z_nkt, c(2,3), FUN=function(x) sum(x,na.rm = TRUE)) %>% as.matrix() %>% t()
    colnames(N.kt) <- c(names.ParetoC, names.BetaC)
#
    # pi gives proportion of mixture components - must add to one across components even when missing some timepoints per mutation
    slices_without_na <- apply(fit$z_nkt, 3, function(slice) sum(apply(slice, 1, function(row) all(!is.na(row))))) # no. of non NA slices for each timepoint
    fit$pi  <- as.list(data.frame(N.kt)/slices_without_na) # Update mixing proportions for each cluster
#
    fit = mobster:::.set_params_Pi(fit, fit$pi)
    
    # PARETO: analytical MLE
    if (fit$fit.tail){
      if ((!log_MMscaled)&(!untf_convolved)){ # normal Pareto tail 
        # fitting scale in non-transformed space was impossible due to numerical problems
        lower = c(1e-4)
        upper = c(15)
              
        for (t in 1:fit$n_t){
          if (i<=genSA_steps){ 
            try( expr={ # if fails, retains previous values
               bf = mobster:::.NLL_Pareto(X$VAF, fit$z_nkt, fit$pi, unlist(fit$scale),t)(unlist(fit$shape)[t])
               MLE.fit = optimization::optim_sa(mobster:::.NLL_Pareto(X$VAF, fit$z_nkt, fit$pi, unlist(fit$scale),t), 
                                               start = c(unlist(fit$shape)[t]),
                                               lower = lower, upper = upper, 
                                               control = list(stopac=3,ac_acc=1e-4))
               
               if ((MLE.fit$function_value - bf) < -1e-4){ # if it actually minimises
                 fit$shape[[1]][t] = MLE.fit$par[1] %>% unname() 
               }
            }) 
          }
          try(expr={
            MLE.fitL = bbmle::mle2(mobster:::.NLL_Pareto(X$VAF, fit$z_nkt, fit$pi, unlist(fit$scale),t),
                                   start = list(shape = unlist(fit$shape)[t]), control = list(trace = 3), skip.hessian = TRUE
                                   );
            
            fit$shape[[1]][t] = as.numeric(unlist(stats4::coef(MLE.fitL))['shape']);  
          })
        }
        fit = mobster:::.set_params_Pareto(fit, fit$shape, fit$scale, fit$times)
      } else if ((!log_MMscaled)&(untf_convolved)){ # convolved Pareto tail 
        # fitting scale in non-transformed space was impossible due to numerical problems
        lower = c(1e-4)
        upper = c(15)
              
        for (t in 1:fit$n_t){
          if (i<=genSA_steps){ 
            try( expr={ # if fails, retains previous values
               bf = mobster:::.NLL_Pareto_convolved(X$VAF, fit$z_nkt, fit$pi, unlist(fit$scale), fit$sigma, t)(unlist(fit$shape)[t])
               MLE.fit = DEoptim::DEoptim(mobster:::.NLL_Pareto_convolved(X$VAF, fit$z_nkt, fit$pi, unlist(fit$scale), fit$sigma, t), 
                                               #start = c(unlist(fit$shape)[t]),
                                               lower = lower, upper = upper, 
                                               control = DEoptim.control(parallelType = 'none',reltol=1e-3,itermax=10*abs(i-0.5) + 5))
                                               #control = list(stopac=3,ac_acc=1e-4))
               
               if ((MLE.fit$optim$bestval - bf) < -1e-4){ # if it actually minimises
                 fit$shape[[1]][t] = MLE.fit$optim$bestmem[1] %>% unname()
               }
            }) 
          }
          try(expr={
            MLE.fitL = bbmle::mle2(mobster:::.NLL_Pareto_convolved(X$VAF, fit$z_nkt, fit$pi, unlist(fit$scale), fit$sigma, t),
                                   start = list(shape = unlist(fit$shape)[t]), control = list(trace = 3), skip.hessian = TRUE
                                   );
            
            fit$shape[[1]][t] = as.numeric(unlist(stats4::coef(MLE.fitL))['shape']);  
          })
        }
        fit = mobster:::.set_params_Pareto(fit, fit$shape, fit$scale, fit$times)
      }
    }
#
    # BETA: numerical MLE or analytical MM
    for (k in 2:fit$K)
    {
        # Compute a functional of the negative logLik, which we minimize
        if ((!log_MMscaled) & untf_convolved){
          for (t in 1:fit$n_t){
            lower = rep(1e-2,2)
            upper = c(75, 350)
            
            if (i<=genSA_steps){ 
              try( expr={ # if fails, retains previous values
                 bf = mobster:::.NLLBetaMixG_convolved(k, X$VAF, fit$z_nkt, fit$pi, fit$sigma, t)(c(fit$a[[k - 1]][t],fit$b[[k - 1]][t]))
                 MLE.fit = DEoptim::DEoptim(mobster:::.NLLBetaMixG_convolved(k, X$VAF, fit$z_nkt, fit$pi, fit$sigma, t), 
                                          lower = lower, upper = upper, control = DEoptim.control(parallelType = 'none',reltol=1e-3,itermax=20*abs(i-0.5) + 10))
#
                 if ((MLE.fit$optim$bestval - min(MLE.fit$member$bestvalit[[1]],bf)) < -1e-4){ # if it actually minimises
                   fit$a[[k - 1]][t] = MLE.fit$optim$bestmem[1] %>% unname() 
                   fit$b[[k - 1]][t] = MLE.fit$optim$bestmem[2]%>% unname() 
               }  
              }) 
            }
            try( expr={ # if fails, retains previous values
              start_list = list(a = fit$a[[k - 1]][t], b = fit$b[[k - 1]][t])
              MLE.fit = bbmle::mle2(
                minuslogl = mobster:::.NLLBetaMixL_convolved(k, X$VAF, fit$z_nkt, fit$pi, fit$sigma, t),
                start = start_list, skip.hessian = TRUE, control=list(trace=3)
                );
              fit$a[[k - 1]][t] = as.numeric(stats4::coef(MLE.fit))[1];
              fit$b[[k - 1]][t] = as.numeric(stats4::coef(MLE.fit))[2];
            }) 
          }
        } else {
          for (t in 1:fit$n_t){
            lower = rep(1e-2,2)
            upper = c(75, 350)
            
            if (i<=genSA_steps){ 
              try( expr={ # if fails, retains previous values
                 bf = mobster:::.NLLBetaMixG(k, X$VAF, fit$z_nkt, fit$pi, t)(c(fit$a[[k - 1]][t],fit$b[[k - 1]][t]))
                 MLE.fit = optimization::optim_sa(mobster:::.NLLBetaMixG(k, X$VAF, fit$z_nkt, fit$pi, t), 
                                                 start= c(fit$a[[k - 1]][t],fit$b[[k - 1]][t]),
                                                 lower = lower, upper = upper, 
                                                 control = list(stopac=3,ac_acc=1e-4))
                 
                 if ((MLE.fit$function_value - bf) < -1e-4){ # if it actually minimises
                   fit$a[[k - 1]][t] = MLE.fit$par[1] %>% unname() 
                   fit$b[[k - 1]][t] = MLE.fit$par[2] %>% unname() 
                 }
              }) 
            }
            try( expr={ # if fails, retains previous values
              start_list = list(a = fit$a[[k - 1]][t], b = fit$b[[k - 1]][t])
              MLE.fit = bbmle::mle2(
                minuslogl = mobster:::.NLLBetaMixL(k, X$VAF, fit$z_nkt, fit$pi, t),
                start = start_list, skip.hessian = TRUE, control=list(trace=3)
                ); 
              fit$a[[k - 1]][t] = as.numeric(stats4::coef(MLE.fit))[1];
              fit$b[[k - 1]][t] = as.numeric(stats4::coef(MLE.fit))[2];
            }) 
          }
        }
      
        fit = mobster:::.set_params_ExpMMbeta(fit, fit$a[k - 1], fit$b[k - 1])
    }
    
    # fitting sigma
    if (log_MMscaled | untf_convolved){
      transformed = ifelse(log_MMscaled,TRUE,FALSE)
      
      if (i<=genSA_steps){ 
        # global minimisation
        print("global")
        
        for (t in 1:fit$n_t){
          try(expr={
            bf = mobster:::.NLL_sigma(X$VAF[[t]],fit,tail,t,transformed=transformed)(fit$sigma[[t]])
            if (Gopt=='DE'){
              MLE.fit = DEoptim::DEoptim(mobster:::.NLL_sigma(X$VAF[[t]],fit,tail,t,transformed=transformed),  # fails here
                                          lower = 1e-3, upper = 0.1, control = DEoptim.control(parallelType = 'none',reltol=1e-3,itermax=20*abs(i-0.5) + 10))
              
              if ((MLE.fit$optim$bestval - min(MLE.fit$member$bestvalit[[1]],bf)) < -1e-4){ # if it actually minimises
                fit$sigma[[t]] = MLE.fit$optim$bestmem %>% unname()
              }
            } else {
             MLE.fit = optimization::optim_sa(mobster:::.NLL_sigma(X$VAF[[t]],fit,tail,t,transformed=transformed), 
                                             start= c(fit$sigma[[t]]),
                                             lower = 1e-3, upper = 0.1, 
                                             control = list(stopac=3,ac_acc=1e-4))
             
             if ((MLE.fit$function_value - bf) < -1e-4){ # if it actually minimises
               fit$sigma[[t]] = MLE.fit$par[1] %>% unname() 
             }
            }
          })
        }
          
        for (t in 1:fit$n_t){
          try( expr={ # if fails, retains previous values
            MLE.fitL = bbmle::mle2(
              minuslogl = mobster:::.NLL_sigma(X$VAF[[t]],fit,tail,t,transformed=transformed),
              start = list(y = fit$sigma[t]), skip.hessian = TRUE, control=list(trace=3)
            );
            fit$sigma[t] = as.numeric(stats4::coef(MLE.fitL));
          })
        }
        
      } else { 
        
        for (t in 1:fit$n_t){
          try( expr={ # if fails, retains previous values
            MLE.fitL = bbmle::mle2(
              minuslogl = mobster:::.NLL_sigma(X$VAF[[t]],fit,tail,t,transformed=transformed),
              start = list(y = fit$sigma[t]), skip.hessian = TRUE, control=list(trace=3)
            );
            fit$sigma[t] = as.numeric(stats4::coef(MLE.fitL));
          })
        } 
      }
    }
        
  return(fit)
  
}
