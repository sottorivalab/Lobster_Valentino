hard_clustered_VAF_CCF <- function(x){
  
  clusters = unique(x$data$cluster)
  non_tail_clusters = clusters[(clusters!='Tail')&(!is.na(clusters))]
  
  hc_attribs = NULL
  
  for (cl in non_tail_clusters){
    
    # VAF standard error on mean 
    vaf_means = x$data %>% filter(cluster==cl) %>% .[['VAF_orig']] %>% lapply(mean, na.rm=TRUE) %>% unlist()
    vaf_sd = x$data %>% filter(cluster==cl) %>% .[['VAF_orig']] %>% lapply(sd, na.rm=TRUE) %>% unlist()
    nobs = x$data %>% filter(cluster==cl) %>% .[['VAF_orig']] %>% lapply(function(x) {sum(!is.na(x))}) %>% unlist()
    vaf_sem = vaf_sd/sqrt(nobs)
    
    hc_attribs =  dplyr::bind_rows(
      hc_attribs,
      tibble::tribble(
        ~cluster, ~type, ~value,
        cl, "VAF_mean", vaf_means, 
        cl, "VAF_sd", vaf_sd,
        cl, "VAF_nobs", nobs,
        cl, "VAF_sem", vaf_sem
      )
    )
  }
  
  if (length(non_tail_clusters)>0){
    hc_attribs = extract_ccf_from_vaf_attribs(hc_attribs,non_tail_clusters,x$n_t)  
    x$hc_attribs = hc_attribs
  } else {
    x$hc_attribs = NULL
  }
  
  return(x)
}

hard_clustered_VAF_CCF_tail <- function(x){
  
  clusters = unique(x$data$cluster)
  clusters = clusters[(!is.na(clusters))]
  
  hc_attribs = NULL
  
  for (cl in clusters){
    
    # VAF standard error on mean 
    vaf_means = x$data %>% filter(cluster==cl) %>% .[['VAF_orig']] %>% lapply(mean, na.rm=TRUE) %>% unlist()
    vaf_sd = x$data %>% filter(cluster==cl) %>% .[['VAF_orig']] %>% lapply(sd, na.rm=TRUE) %>% unlist()
    nobs = x$data %>% filter(cluster==cl) %>% .[['VAF_orig']] %>% lapply(function(x) {sum(!is.na(x))}) %>% unlist()
    vaf_sem = vaf_sd/sqrt(nobs)
    
    hc_attribs =  dplyr::bind_rows(
      hc_attribs,
      tibble::tribble(
        ~cluster, ~type, ~value,
        cl, "VAF_mean", vaf_means, 
        cl, "VAF_sd", vaf_sd,
        cl, "VAF_nobs", nobs,
        cl, "VAF_sem", vaf_sem
      )
    )
  }
  
  if (length(clusters)>0){
    hc_attribs = extract_ccf_from_vaf_attribs(hc_attribs,clusters,x$n_t)  
    x$hc_attribs = hc_attribs
  } else {
    x$hc_attribs = NULL
  }
  
  return(x)
}


densities_VAF_CCF <- function(x){
  
  clusters = unique(x$data$cluster)
  non_tail_clusters = clusters[clusters!='Tail']
  
  dens_attribs = NULL
  
  for (cl in non_tail_clusters){
    
    # VAF standard error on mean 
    vaf_means = x$Clusters %>% filter(cluster==cl,type=='Mean') %>% select(fit.value) %>% unlist() %>% unname()
    vaf_sd = x$Clusters %>% filter(cluster==cl,type=='Variance') %>% select(fit.value) %>% unlist() %>% sapply(FUN= function(x) {sqrt(x)}) %>% unname()
    nobs = x$data %>% filter(cluster==cl) %>% .[['VAF_orig']] %>% lapply(function(x) {sum(!is.na(x))}) %>% unlist() %>% unname()
    vaf_sem = vaf_sd/sqrt(nobs) %>% unname()
    
    dens_attribs =  dplyr::bind_rows(
      dens_attribs,
      tibble::tribble(
        ~cluster, ~type, ~value,
        cl, "VAF_mean", vaf_means, 
        cl, "VAF_sd", vaf_sd,
        cl, "VAF_nobs", nobs,
        cl, "VAF_sem", vaf_sem
      )
    )
  }
  
  if (length(non_tail_clusters)>0){
    dens_attribs = extract_ccf_from_vaf_attribs(dens_attribs,non_tail_clusters,x$n_t)  
    x$dens_attribs = dens_attribs
  } else {
    x$dens_attribs = NULL
  }
  
  return(x)
}

extract_ccf_from_vaf_attribs <- function(attribs,non_tail_clusters,n_t){

  if ("C1" %in% non_tail_clusters){
    vaf_mean_c1 = filter(attribs, cluster=="C1", type=="VAF_mean") %>% select(value) %>% unlist()  
    vaf_sem_c1 = filter(attribs, cluster=="C1", type=="VAF_sem") %>% select(value) %>% unlist()
    
    for (cl in non_tail_clusters){
      # C1 should be the cluster with the highest mean
      # assuming all clusters have same ploidy
      if (cl != 'C1'){
        vaf_sem = filter(attribs, cluster==cl, type=="VAF_sem") %>% select(value) %>% unlist()
        vaf_means = filter(attribs, cluster==cl, type=="VAF_mean") %>% select(value) %>% unlist()
        
        ccf = ((vaf_means)/(vaf_mean_c1)) %>% unname()
        ccf_err = (ccf*sqrt((vaf_sem/vaf_means)^2 + (vaf_sem_c1/vaf_mean_c1)^2))  %>% unname() # standard error propagation
        
      } else { 
        ccf = rep(1.0,n_t)
        ccf_err = rep(NA,n_t)
      }
    
      attribs =  dplyr::bind_rows(
        attribs,
        tibble::tribble(
          ~cluster, ~type, ~value,
          cl, "CCF", ccf, 
          cl, "CCF_err", ccf_err
        )
      )
    }
  }
  
  return(attribs)
}  

extract_ccf_from_vafs_tibble <- function(unscaled_vaf,fit){

  vaf_mean_c1 = filter(fit$hc_attribs, cluster=="C1", type=="VAF_mean") %>% select(value) %>% unlist()  
  vaf_sem_c1 = filter(fit$hc_attribs, cluster=="C1", type=="VAF_sem") %>% select(value) %>% unlist()
  
  fit$data$CCF = fit$data$VAF_orig %>% mutate_all(~ replace(., TRUE, NA)) # initialise copy of VAF tibble with NA 
  fit$data$CCF_errs = fit$data$VAF_orig %>% mutate_all(~ replace(., TRUE, NA)) # initialise copy of VAF tibble with NA 
  for (t in 1:fit$n_t){
    fit$data$CCF[[t]] = unscaled_vaf[[t]]/(vaf_mean_c1[t])
    fit$data$CCF_errs[[t]] = fit$data$CCF[[t]]*((vaf_sem_c1[t])/(vaf_mean_c1[t])) # standard error propagation
  }
      
  return(fit)
}  

ccf_to_vaf <- function(unscaled_vaf,fit,t){

  vaf_mean_c1 = filter(fit$hc_attribs, cluster=="C1", type=="VAF_mean") %>% select(value) %>% unlist()  
  ccf = unscaled_vaf/(vaf_mean_c1[t])
      
  return(ccf)
}  

