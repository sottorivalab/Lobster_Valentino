.Mstep <- function(i, X, fit, MM_scale_a, MM_scale_b, 
                  genSA_steps, log_MMscaled, untf_convolved, Gopt)
{
  transformed = log_MMscaled
  
  if (log_MMscaled){ # expMM 
     minimisation('beta','expMM', X$VAF, fit, transformed, Gopt, i, genSA_steps)
     if (fit$fit.tail){
       # optimising scale with fixed shape
       minimisation('pareto','expMM_scale', X$VAF, fit, transformed, Gopt, i, genSA_steps)
       # optimising shape with fixed scale
  browser()
       minimisation('pareto','expMM_shape', X$VAF, fit, transformed, Gopt, i, genSA_steps)
     }   
  } else if ((!log_MMscaled) & untf_convolved){ # convolved
     minimisation('beta','convolved', X$VAF, fit, transformed, Gopt, i, genSA_steps)
     if (fit$fit.tail){
       minimisation('pareto','convolved', X$VAF, fit, transformed, Gopt, i, genSA_steps)
     }
  } else { # normal Pareto tail 
     minimisation('beta','simple', X$VAF, fit, transformed, Gopt, i, genSA_steps)
     if (fit$fit.tail){
       minimisation('pareto','simple', X$VAF, fit, transformed, Gopt, i, genSA_steps)
     }
  }
  
  # fitting sigma
  if (log_MMscaled | untf_convolved){
    minimisation('sigma','', X$VAF, fit, transformed, Gopt, i, genSA_steps)
  }  
  
  return(fit)  
}

fit_beta_simpleL <- function(vaf, fit, transformed) {

  for (k in 2:fit$K)
  {
    for (t in 1:fit$n_t){
       try( expr={ # if fails, retains previous values
         start_list = list(a = fit$a[[k - 1]][t], b = fit$b[[k - 1]][t])
         MLE.fit = bbmle::mle2(
           minuslogl = mobster:::.NLLBetaMixL(k, vaf, fit$z_nkt, fit$pi, t),
           start = start_list, skip.hessian = TRUE, control=list(trace=3)
           );
         fit$a[[k - 1]][t] = as.numeric(stats4::coef(MLE.fit))[1];
         fit$b[[k - 1]][t] = as.numeric(stats4::coef(MLE.fit))[2];
       }) 
    }
    fit = mobster:::.set_params_ExpMMbeta(fit, fit$a[k - 1], fit$b[k - 1])
  }

  return(fit)
}  

fit_beta_convolvedL <- function(vaf, fit, transformed) {
  
  for (k in 2:fit$K)
  {
     for (t in 1:fit$n_t){
      try( expr={ # if fails, retains previous values
        start_list = list(a = fit$a[[k - 1]][t], b = fit$b[[k - 1]][t])
        MLE.fit = bbmle::mle2(
          minuslogl = mobster:::.NLLBetaMixL_convolved(k, vaf, fit$z_nkt, fit$pi, fit$sigma, t),
          start = start_list, skip.hessian = TRUE, control=list(trace=3)
          );
        fit$a[[k - 1]][t] = as.numeric(stats4::coef(MLE.fit))[1];
        fit$b[[k - 1]][t] = as.numeric(stats4::coef(MLE.fit))[2];
      }) 
    }
    fit = mobster:::.set_params_ExpMMbeta(fit, fit$a[k - 1], fit$b[k - 1])
  }

  return(fit)
}  

fit_beta_expMML <- function(vaf, fit, transformed) {
  
  for (k in 2:fit$K)
  {
     for (t in 1:fit$n_t){
       start_list = list(a = fit$a[[k - 1]][t], b = fit$b[[k - 1]][t])
       try( expr={ # if fails, retains previous values 
         MLE.fit = bbmle::mle2(
           minuslogl = mobster:::.NLL_ExpMMBetaMixL_convolved(k, vaf[[t]], fit$z_nkt, fit$pi, fit$MM_scale$MM_scale_a, fit$MM_scale$MM_scale_b, fit$sigma, t),
           start = start_list, skip.hessian = TRUE, control=list(trace=3));
         fit$a[[k - 1]][t] = as.numeric(stats4::coef(MLE.fit))[1];
         fit$b[[k - 1]][t] = as.numeric(stats4::coef(MLE.fit))[2];
       })
     }
     fit = mobster:::.set_params_ExpMMbeta(fit, fit$a[k - 1], fit$b[k - 1])
  }

  return(fit)
}  

fit_pareto_simpleL <- function(vaf, fit, transformed) {
  
  for (t in 1:fit$n_t){
    try(expr={
      MLE.fitL = bbmle::mle2(mobster:::.NLL_Pareto(vaf, fit$z_nkt, fit$pi, unlist(fit$scale),t),
                             start = list(shape = unlist(fit$shape)[t]), control = list(trace = 3), skip.hessian = TRUE
                             );
      
      fit$shape[[1]][t] = as.numeric(unlist(stats4::coef(MLE.fitL))['shape']);  
    })
  }

  fit = mobster:::.set_params_Pareto(fit, fit$shape, fit$scale, fit$times)
  
  return(fit)
}  

fit_pareto_convolvedL <- function(vaf, fit, transformed) {
  
  for (t in 1:fit$n_t){
    try(expr={
      MLE.fitL = bbmle::mle2(mobster:::.NLL_Pareto_convolved(vaf, fit$z_nkt, fit$pi, unlist(fit$scale), fit$sigma, t),
                             start = list(shape = unlist(fit$shape)[t]), control = list(trace = 3), skip.hessian = TRUE
                             );
      
      fit$shape[[1]][t] = as.numeric(unlist(stats4::coef(MLE.fitL))['shape']);  
    })
  }

  fit = mobster:::.set_params_Pareto(fit, fit$shape, fit$scale, fit$times)
  
  return(fit)
}  

fit_pareto_expMM_shapeL <- function(vaf, fit, transformed) {

  for (t in 1:fit$n_t){
    try(expr={
       for (t in 1:fit$n_t){
            MLE.fitL = bbmle::mle2(
              minuslogl = mobster:::.NLL_ExpMMtailL_shape(vaf[[t]], fit$z_nkt, fit$pi, fit$MM_scale$MM_scale_a, fit$MM_scale$MM_scale_b, fit$sigma, unlist(fit$scale), t),
              start = list(shape = unlist(fit$shape)[t]), 
              control = list(trace = 3,maxit=100, REPORT=1), skip.hessian = TRUE
            )
          fit$shape[[1]][t] = as.numeric(unlist(stats4::coef(MLE.fitL))['shape']);
       }
    })
  }

  fit = mobster:::.set_params_ExpMMtail(fit, fit$shape, fit$scale) 
  
  return(fit)
}  

fit_pareto_expMM_scaleL <- function(vaf, fit, transformed) {
  
  try(expr={
      MLE.fitL = bbmle::mle2(
        minuslogl = mobster:::.NLL_ExpMMtailL_scale(vaf, fit$z_nkt, fit$pi, fit$times, fit$MM_scale$MM_scale_a, fit$MM_scale$MM_scale_b, fit$sigma, unlist(fit$shape)),
        start = list(scale = unlist(fit$scale)), 
        control = list(trace = 3,maxit=100, REPORT=1), skip.hessian = TRUE
      )
    fit$scale[[1]] = as.numeric(stats4::coef(MLE.fitL)['scale']);
  })
  
  fit = mobster:::.set_params_ExpMMtail(fit, fit$shape, fit$scale) 

  return(fit)
}  

fit_sigmaL <- function(vaf, fit, transformed) {
  
  for (t in 1:fit$n_t){
    try( expr={ # if fails, retains previous values
      MLE.fitL = bbmle::mle2(
        minuslogl = mobster:::.NLL_sigma(vaf[[t]],fit,fit$fit.tail,t,transformed=transformed),
        start = list(y = fit$sigma[t]), skip.hessian = TRUE, control=list(trace=3)
      );
      fit$sigma[t] = as.numeric(stats4::coef(MLE.fitL));
    })
  } 
  
  return(fit)
}

fit_sigmaG <- function(vaf, fit, transformed, Gopt, i) {
  
  for (t in 1:fit$n_t){
    bf = mobster:::.NLL_sigma(vaf[[t]],fit,fit$fit.tail,t,transformed=transformed)(fit$sigma[[t]])
    if (Gopt=='DE'){
      MLE.fit = DEoptim::DEoptim(mobster:::.NLL_sigma(vaf[[t]],fit,fit$fit.tail,t,transformed=transformed), 
                                  lower = 1e-3, upper = 0.1, control = DEoptim.control(parallelType = 'none',reltol=1e-3,itermax=10*i))
      
      if ((MLE.fit$optim$bestval - min(MLE.fit$member$bestvalit[[1]],bf)) < -1e-4){ # if it actually minimises
        fit$sigma[[t]] = MLE.fit$optim$bestmem %>% unname()
      }
     } else {
     MLE.fit = optimization::optim_sa(mobster:::.NLL_sigma(vaf[[t]],fit,fit$fit.tail,t,transformed=transformed), 
                                     start= c(fit$sigma[[t]]),
                                     lower = 1e-3, upper = 0.1, 
                                     control = list(stopac=3,ac_acc=1e-4))
     
     if ((MLE.fit$function_value - bf) < -1e-4){ # if it actually minimises
       fit$sigma[[t]] = MLE.fit$par[1] %>% unname() 
     }
   }
  }
  
  return(fit)
}

fit_beta_simpleG <- function(vaf, fit, transformed, Gopt, i) {
  
  lower = rep(1e-2,2)
  upper = c(75, 350)
  
  for (k in 2:fit$K)
  {
     for (t in 1:fit$n_t){
       if (i<=genSA_steps){ 
         try( expr={ # if fails, retains previous values
            bf = mobster:::.NLLBetaMixG(k, vaf, fit$z_nkt, fit$pi, t)(c(fit$a[[k - 1]][t],fit$b[[k - 1]][t]))
            MLE.fit = optimization::optim_sa(mobster:::.NLLBetaMixG(k, vaf, fit$z_nkt, fit$pi, t), 
                                            start= c(fit$a[[k - 1]][t],fit$b[[k - 1]][t]),
                                            lower = lower, upper = upper, 
                                            control = list(stopac=3,ac_acc=1e-4))
            
            if ((MLE.fit$function_value - bf) < -1e-4){ # if it actually minimises
              fit$a[[k - 1]][t] = MLE.fit$par[1] %>% unname() 
              fit$b[[k - 1]][t] = MLE.fit$par[2] %>% unname() 
            }
         }) 
       }
     }
     
     fit = mobster:::.set_params_ExpMMbeta(fit, fit$a[k - 1], fit$b[k - 1])
  }

  return(fit)
}  

fit_beta_convolvedG <- function(vaf, fit, transformed, Gopt, i) {
  
  lower = rep(1e-2,2)
  upper = c(75, 350)
    
  for (k in 2:fit$K)
  {
     for (t in 1:fit$n_t){
       if (i<=genSA_steps){ 
         try( expr={ # if fails, retains previous values
            bf = mobster:::.NLLBetaMixG_convolved(k, vaf, fit$z_nkt, fit$pi, fit$sigma, t)(c(fit$a[[k - 1]][t],fit$b[[k - 1]][t]))
            MLE.fit = DEoptim::DEoptim(mobster:::.NLLBetaMixG_convolved(k, vaf, fit$z_nkt, fit$pi, fit$sigma, t), 
                                     lower = lower, upper = upper, control = DEoptim.control(parallelType = 'none',reltol=1e-3,itermax=20*i))
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
           minuslogl = mobster:::.NLLBetaMixL_convolved(k, vaf, fit$z_nkt, fit$pi, fit$sigma, t),
           start = start_list, skip.hessian = TRUE, control=list(trace=3)
           );
         fit$a[[k - 1]][t] = as.numeric(stats4::coef(MLE.fit))[1];
         fit$b[[k - 1]][t] = as.numeric(stats4::coef(MLE.fit))[2];
       }) 
     }

     fit = mobster:::.set_params_ExpMMbeta(fit, fit$a[k - 1], fit$b[k - 1])
  }

  return(fit)
}  

fit_beta_expMMG <- function(vaf, fit, transformed, Gopt, i) {
  
  lower = rep(1e-2,2)
  upper = c(75, 350)
  
  for (k in 2:fit$K)
  {
     for (t in 1:fit$n_t){
       bf = mobster:::.NLL_ExpMMBetaMixG_convolved(k, vaf[[t]], fit$z_nkt, fit$pi, 
                                                   fit$MM_scale$MM_scale_a, fit$MM_scale$MM_scale_b, fit$sigma, 
                                                   t)(c(fit$a[[k - 1]][t],fit$b[[k - 1]][t]))
      if (Gopt=='DE'){
       MLE.fit = DEoptim::DEoptim(mobster:::.NLL_ExpMMBetaMixG_convolved(k, vaf[[t]], fit$z_nkt, fit$pi, fit$MM_scale$MM_scale_a, fit$MM_scale$MM_scale_b, fit$sigma, t), 
                                  lower = lower, upper = upper, control = DEoptim.control(parallelType = 'none',reltol=1e-3,itermax=20*i))
   #
       if ((MLE.fit$optim$bestval - min(MLE.fit$member$bestvalit[[1]],bf)) < -1e-4){ # if it actually minimises
         fit$a[[k - 1]][t] = MLE.fit$optim$bestmem[1] %>% unname() 
         fit$b[[k - 1]][t] = MLE.fit$optim$bestmem[2]%>% unname() 
       }
     }
   }

   fit = mobster:::.set_params_ExpMMbeta(fit, fit$a[k - 1], fit$b[k - 1])
  }

  return(fit)
}  

fit_pareto_simpleG <- function(vaf, fit, transformed, Gopt, i) {
  
  lower = c(1e-4)
  upper = c(15)
        
  for (t in 1:fit$n_t){
    try( expr={ # if fails, retains previous values
       bf = mobster:::.NLL_Pareto(vaf, fit$z_nkt, fit$pi, unlist(fit$scale),t)(unlist(fit$shape)[t])
       MLE.fit = optimization::optim_sa(mobster:::.NLL_Pareto(vaf, fit$z_nkt, fit$pi, unlist(fit$scale),t), 
                                       start = c(unlist(fit$shape)[t]),
                                       lower = lower, upper = upper, 
                                       control = list(stopac=3,ac_acc=1e-4))
       
       if ((MLE.fit$function_value - bf) < -1e-4){ # if it actually minimises
         fit$shape[[1]][t] = MLE.fit$par[1] %>% unname() 
       }
    }) 
  }  
      
  fit = mobster:::.set_params_Pareto(fit, fit$shape, fit$scale, fit$times)
  
  return(fit)
}  

fit_pareto_convolvedG <- function(vaf, fit, transformed, Gopt, i) {
  
  lower = c(1e-4)
  upper = c(15)
        
  for (t in 1:fit$n_t){
    try( expr={ # if fails, retains previous values
       bf = mobster:::.NLL_Pareto_convolved(vaf, fit$z_nkt, fit$pi, unlist(fit$scale), fit$sigma, t)(unlist(fit$shape)[t])
       MLE.fit = DEoptim::DEoptim(mobster:::.NLL_Pareto_convolved(vaf, fit$z_nkt, fit$pi, unlist(fit$scale), fit$sigma, t), 
                                       #start = c(unlist(fit$shape)[t]),
                                       lower = lower, upper = upper, 
                                       control = DEoptim.control(parallelType = 'none',reltol=1e-3,itermax=10*i))
                                       #control = list(stopac=3,ac_acc=1e-4))
       
       if ((MLE.fit$optim$bestval - bf) < -1e-4){ # if it actually minimises
         fit$shape[[1]][t] = MLE.fit$optim$bestmem[1] %>% unname()
       }
    }) 
  }
  
  fit = mobster:::.set_params_Pareto(fit, fit$shape, fit$scale, fit$times)
  
  return(fit)
}  

fit_pareto_expMM_shapeG <- function(vaf, fit, transformed, Gopt, i) {
  
  lower = 1e-4
  upper = 15
  
 # optimising shape with fixed scale
  if (Gopt=='DE'){
   for (t in 1:fit$n_t){
      bf = (mobster:::.NLL_ExpMMtailG_shape(vaf[[t]], fit$z_nkt, fit$pi, fit$MM_scale$MM_scale_a, fit$MM_scale$MM_scale_b, fit$sigma, unlist(fit$scale), t))(unlist(fit$shape)[t])
      MLE.fit = DEoptim::DEoptim(mobster:::.NLL_ExpMMtailG_shape(vaf[[t]], fit$z_nkt, fit$pi, fit$MM_scale$MM_scale_a, fit$MM_scale$MM_scale_b, fit$sigma, unlist(fit$scale), t), 
                                  lower = lower, upper = upper, control = DEoptim.control(parallelType = 'none',reltol=1e-3,itermax=10*i))

      if ((MLE.fit$optim$bestval - min(MLE.fit$member$bestvalit[[1]],bf)) < -1e-4){ # if it actually minimises
        fit$shape[[1]][t] = MLE.fit$optim$bestmem[1] %>% unname() 
      }
   }
  } else {
   bf = (mobster:::.NLL_ExpMMtailG_shape(vaf[[t]], fit$z_nkt, fit$pi, fit$MM_scale$MM_scale_a, fit$MM_scale$MM_scale_b, fit$sigma, unlist(fit$scale), t))(unlist(fit$shape)[t])
   MLE.fit = optimization::optim_sa(mobster:::.NLL_ExpMMtailG_shape(vaf[[t]], fit$z_nkt, fit$pi, 
                                                                    fit$MM_scale$MM_scale_a, fit$MM_scale$MM_scale_b, fit$sigma, unlist(fit$scale), t), 
                                   start= c(unlist(fit$shape)[t]),
                                   lower = lower, upper = upper, 
                                   control = list(stopac=3,ac_acc=1e-4))
   
   if ((MLE.fit$function_value - bf) < -1e-4){ # if it actually minimises
     fit$shape[[1]][t] = MLE.fit$par[1] %>% unname() 
   }
  }
  
  fit = mobster:::.set_params_ExpMMtail(fit, fit$shape, fit$scale) 
  
  return(fit)
}  

fit_pareto_expMM_scaleG <- function(vaf, fit, transformed, Gopt, i) {
  
  lower = 1e-4
  upper = 0.5
  
 # optimising scale with fixed shape
 bf = (mobster:::.NLL_ExpMMtailG_scale(vaf, fit$z_nkt, fit$pi, fit$times, fit$MM_scale$MM_scale_a, fit$MM_scale$MM_scale_b, fit$sigma, unlist(fit$shape)))(unlist(fit$scale))
 
  if (Gopt=='DE'){
   MLE.fit = DEoptim::DEoptim(mobster:::.NLL_ExpMMtailG_scale(vaf, fit$z_nkt, fit$pi, fit$times, fit$MM_scale$MM_scale_a, fit$MM_scale$MM_scale_b, fit$sigma, unlist(fit$shape)), 
                               lower = lower, upper = upper, control = DEoptim.control(parallelType = 'none',reltol=1e-3,itermax=10*i))

   if ((MLE.fit$optim$bestval - min(MLE.fit$member$bestvalit[[1]],bf)) < -1e-4){ # if it actually minimises
     fit$scale[[1]] = MLE.fit$optim$bestmem[1] %>% unname() 
   }
  } else {
   MLE.fit = optimization::optim_sa(mobster:::.NLL_ExpMMtailG_scale(vaf, fit$z_nkt, fit$pi, fit$times, 
                                                                    fit$MM_scale$MM_scale_a, fit$MM_scale$MM_scale_b, fit$sigma, unlist(fit$shape)), 
                                   start= c(unlist(fit$scale)),
                                   lower = lower, upper = upper, 
                                   control = list(stopac=3,ac_acc=1e-4))
   
   if ((MLE.fit$function_value - bf) < -1e-4){ # if it actually minimises
     fit$scale[[1]] = MLE.fit$par[1] %>% unname() 
   }
  }
  
  fit = mobster:::.set_params_ExpMMtail(fit, fit$shape, fit$scale) 
  
  return(fit)
}  

global_min_func <- function(comp_str, sit_str){
  
  if (comp_str == 'sigma'){
    
    func = fit_sigmaG
    
  } 
  else if (comp_str == 'pareto'){
    
      if (sit_str == 'simple'){
        func = fit_pareto_simpleG
        
      } else if (sit_str == 'convolved'){
        func = fit_pareto_convolvedG
        
      } else if (sit_str == 'expMM_scale'){
        func = fit_pareto_expMM_scaleG
        
      } else if (sit_str == 'expMM_shape'){
        func = fit_pareto_expMM_shapeG
        
      } else {simpleError('Invalid sit_str')}
    
  } else if (comp_str == 'beta'){
    
      if (sit_str == 'simple'){
        func = fit_beta_simpleG
        
      } else if (sit_str == 'convolved'){
        func = fit_beta_convolvedG
        
      } else if (sit_str == 'expMM'){
        func = fit_beta_expMMG
        
      } else {simpleError('Invalid sit_str')}
    
  } else {simpleError('Invalid comp_str')}
  
  return(func)
}

local_min_func <- function(comp_str, sit_str){
  
  if (comp_str == 'sigma'){
    
    func = fit_sigmaL
    
  } 
  else if (comp_str == 'pareto'){
    
      if (sit_str == 'simple'){
        func = fit_pareto_simpleL
        
      } else if (sit_str == 'convolved'){
        func = fit_pareto_convolvedL
        
      } else if (sit_str == 'expMM_scale'){
        func = fit_pareto_expMM_scaleL
        
      } else if (sit_str == 'expMM_shape'){
        func = fit_pareto_expMM_shapeL
        
      } else {simpleError('Invalid sit_str')}
    
  } else if (comp_str == 'beta'){
    
      if (sit_str == 'simple'){
        func = fit_beta_simpleL
        
      } else if (sit_str == 'convolved'){
        func = fit_beta_convolvedL
        
      } else if (sit_str == 'expMM'){
        func = fit_beta_expMML
        
      } else {simpleError('Invalid sit_str')}
    
  } else {simpleError('Invalid comp_str')}
  
  return(func)
}

minimisation <- function(comp_str, sit_str, vaf, fit, transformed, Gopt, i, genSA_steps){
  
  if (i<=genSA_steps){ 
    gfunc = global_min_func(comp_str,sit_str)
    fit = gfunc(vaf, fit, transformed, Gopt, i)  
  } 
  
  lfunc = local_min_func(comp_str,sit_str)
  fit = lfunc(vaf, fit, transformed)  
  
  return(fit)
}

