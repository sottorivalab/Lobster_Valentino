#' Filter MOBSTER output clusters.
#' 
#' @description This function can filter out the clusters computed by MOBSTER
#' based on two criteria: the mixing proportion value,  the number of mutations
#' assigned and the variance of the Beta clusters. 
#' 
#' For all criteria a scalar should be given as input. The return object
#' will contain only the clusters that pass all filters. If any cluster is dropped
#' the latent variables are re-computed, as well as the clustering assignments and the
#' mixing proportions (all mutations will be still assigned after clusters' removal).
#' 
#' @param x A MOBSTER fit object.
#' @param pi_cutoff The cutoff on the mixing proportions, default is 0.02.
#' @param N_cutoff The cutoff on the number of mutations assigned to a cluster, default is 10.  
#' @param Beta_variance_cutoff Minimum variance for a Beta peak.
#' @param verbose If outputs should be reported to screen or not, default is no.
#'
#' @return A MOBSTER fit object where clusters are larger than \code{pi_cutoff} and contain
#' at least \code{N_cutoff}. If no such cluster exists an error is generated.
#' 
#' @export
#' 
#' @examples
#' data('fit_example', package = 'mobster')
#'
#' # Does not change anything (no filter triggered)
#' choose_clusters(fit_example$best)
#' 
#' # Remove one Beta component because it has less than 100 points (renders the fit very poor)
#' choose_clusters(fit_example$best, N_cutoff = 100)
choose_clusters = function(x, 
                           pi_cutoff = 0.02,
                           N_cutoff = 10,
                           Beta_variance_cutoff = 0.00001,
                           log_MMscaled = FALSE,
                           untf_convolved = FALSE,
                           MM_scale_a = NULL,
                           MM_scale_b = NULL,
                           verbose = FALSE)
{
  mobster:::is_mobster_fit(x)
  
  if(verbose)
  {
    cli::cli_rule("Filtering MOBSTER clusters")
    cli::cli_text(mobster:::m_inf(paste0("pi > {.field {", pi_cutoff, '}} (mixing proportions)')))
    cli::cli_text(mobster:::m_inf(paste0(" n > {.field {", N_cutoff, '}} (mutations per cluster)')))
    cli::cli_text(mobster:::m_inf(paste0(" v > {.field {", Beta_variance_cutoff, '}} (minimum Beta variance)')))
  }
  
  # Cluster size - remove clusters smaller than pi_cutoff, or less than N_cutoff mutations
  pass_clusters_Ncutoff = x$N.k > N_cutoff
  
  # Remove beta clusters with variance too small
  # requirement should be sum of variances across timepoints weighted by pi_t/n_t 
  pi = mobster:::.params_Pi(x)
  if (is.null(MM_scale_b)){varb = x$Clusters %>% dplyr::filter(type == 'Variance', cluster != 'Tail')}
  pass_clusters_Bvcutoff = list()
  if (is.null(MM_scale_b)){
    for (i in varb$cluster){
      pass_clusters_Bvcutoff[i] = (((varb %>% dplyr::filter(cluster==i) %>% select(fit.value) %>% unlist() %>% as.vector()) %*%
                                      as.vector(x$pi[[i]])) / x$n_t) > Beta_variance_cutoff  
    }
  }  
  
  tab_clusters = tibble::tibble(cluster = names(pi), N = x$N.k[names(pi)])
  tab_clusters$F1 = pass_clusters_Ncutoff
  if (is.null(MM_scale_b)){ 
    tab_clusters = tab_clusters %>% 
      dplyr::left_join(
        data.frame(stringsAsFactors = FALSE, 
                   cluster = names(pass_clusters_Bvcutoff),  #Beta_var = varb$fit.value
                   F3 = unlist(pass_clusters_Bvcutoff)),
        by = 'cluster'
        ) %>%
      dplyr::select(cluster, N, starts_with('F')) %>% replace(., is.na(.), TRUE)
    tab_clusters$retain = apply(tab_clusters, 1, function(x) all(as.logical(x[3:length(x)]))) 
  }
  else{
    tab_clusters = tab_clusters %>% dplyr::select(cluster, N, starts_with('F')) %>% replace(., is.na(.), TRUE)
    tab_clusters$retain = apply(tab_clusters, 1, function(x) all(as.logical(x[3:length(x)])))
  }
  if (!(x$fit.tail)){tab_clusters[tab_clusters$cluster=='Tail',]$retain = FALSE}
  
  # CLusters that will be cancelled
  clusters_to_cancel = tab_clusters$cluster[which(!tab_clusters$retain)]
    
  # Report to screen
  if(verbose) print(tab_clusters %>% arrange(desc(retain)))
  
  # Report easy cases
  if(sum(!tab_clusters$retain) == nrow(tab_clusters)) stop("All clusters filtered out, this seems an error.")
  if(sum(tab_clusters$retain) == 0) return(x)
  
  # Clusters that remain, which we remove from a copy of x
  remaining_clusters = tab_clusters$cluster[which(tab_clusters$retain)]
  remaining_beta_clusters = remaining_clusters[grepl('C', remaining_clusters)]
  
  # We modify a copy of the input
  y = x
  
  # Remove from the clusters table any information regarding old clusters,
  # preserving the tail's proportion as a valid entry
  y$Clusters = y$Clusters %>% dplyr::filter(cluster %in% remaining_clusters)
  y$fit.tail = ((x$fit.tail)&('Tail' %in% remaining_clusters))
  
  if(!(y$fit.tail)) 
    y$Clusters = y$Clusters %>%
    dplyr::bind_rows(
      tibble::tribble(
        ~cluster, ~type, ~fit.value, ~init.value,
        'Tail', "Mixing proportion", as.list(rep(0,x$n_t)), rep(NA,x$n_t)
      )
    )
  
  # Update counts for Beta components
  y$Kbeta = length(remaining_beta_clusters)
  y$K = y$Kbeta + 1
  
  # Renormalize the mixing proportions -- adding 0% Tail if it's no longer fit
  # and ensure that the ordering is correct 
  y$pi = x$pi[remaining_clusters]
  if(!(y$fit.tail)){
    y$pi[['Tail']] = rep(0,x$n_t)
  }  
  norms = Reduce(`+`, y$pi)
  for (i in 1:length(y$pi))
    for (j in 1:length(y$pi[[i]])){
      y$pi[[i]][j] = y$pi[[i]][j]/norms[j]  
    }
    
  if(!(y$fit.tail)) {
    y$pi = y$pi[c('Tail', remaining_clusters)]
  }
  
  y = mobster:::.set_params_Pi(y, y$pi)
  
  # If we need to cancel a tail
  if(!(y$fit.tail)) 
  {
    y$scale = list(NA)
    y$shape = list(rep(NA,x$n_t))
    # This should be useless
    if (!log_MMscaled){
      y = mobster:::.set_params_Pareto(y, y$scale, y$shape, y$times)
    } else {
      y = mobster:::.set_params_ExpMMtail(y, y$scale, y$shape)
    }
  }
  
  # For the Beta that remains, we need to subset their parameters
  if(length(remaining_beta_clusters) > 0)
  {
    y$a = y$a[remaining_beta_clusters]
    y$b = y$b[remaining_beta_clusters]
    if (!log_MMscaled){
      y = mobster:::.set_params_Beta(y, y$a, y$b, y$times)
    } else {
      y = mobster:::.set_params_ExpMMbeta(y, y$a, y$b)
    }
  }
  
  new_lv = mobster:::latent_vars(y, 
                  log_MMscaled = log_MMscaled,
                  MM_scale_a = MM_scale_a,
                  MM_scale_b = MM_scale_b,
                  untf = untf_convolved
                  )
  
  y$NLL = new_lv$NLL
  y$z_nkt = new_lv$z_nkt
  y$pdf.w = new_lv$pdf.w
  
  # Update clustering assignments and summary numbers
  y$data$cluster = mobster:::latent_vars_hard_assignments(new_lv)
  
  y$N.k = rep(0, y$K)
  names(y$N.k) = names(y$pi)
  obFreq = table(y$data$cluster)
  
  y$N.k[names(obFreq)] = obFreq
  
  # Update scores for model selection   
  y$scores = mobster:::latent_vars_scores(
    mobster:::latent_vars(y,
                  log_MMscaled = log_MMscaled,
                  MM_scale_a = MM_scale_a, 
                  MM_scale_b = MM_scale_b,
                  untf = untf_convolved), # Extract latent variables
    y$K,
    y$n_t,
    y$fit.tail,
    y$data$cluster,
    (log_MMscaled|untf_convolved))
  
  if(verbose)
  {
    pio::pioTit("Selected clusters")
    print(y)
  }
  
  y
}

get_cluster_mode = function(x, MM_scale_a, MM_scale_b)
{
  domain = seq(0, 1, 0.001)
  densities = suppressWarnings(mobster:::template_density(
    x,
    x.axis = domain,
    binwidth = 0.001,
    reduce = TRUE,
    MM_scale_a = MM_scale_a,
    MM_scale_b = MM_scale_b,
    log_MMscaled = TRUE
  ))
  mode_x = list()
  for (i in unique(densities$cluster)){
    mode_idx = densities %>% dplyr::filter(cluster==i) %>% select(y) %>%
      unlist() %>% which.max()
    mode_x[[i]] <- domain[mode_idx]
  }
  
  return(mode_x)
}

get_cluster_mean_and_var = function(x,log_MMscaled)
{
  # calculating on untransformed template density - these are kept untransformed
  binwidth = 0.0001
  domain = seq(0, 1, binwidth)
  
  densities = suppressWarnings(mobster:::template_density(
    x,
    x.axis = domain,
    binwidth = binwidth,
    reduce = TRUE,
    MM_scale_a = MM_scale_a,
    MM_scale_b = MM_scale_b,
    log_MMscaled = log_MMscaled,
    untf = TRUE
  ))
  
  # create times column and rbind densities together
  densities1 = lapply(seq_along(densities),function(y,n,i){ mutate(y[[i]],
                times = as.integer(substr(n[[i]],nchar(n[[i]]),nchar(n[[i]])))) },y=densities,n=names(densities))
  densities <- do.call(rbind, densities1)
  
  for (i in unique(densities$cluster)){
    clusteri = densities %>% filter(cluster==i)
    mean_x = rep(NA,x$n_t)
    var_x = rep(NA,x$n_t)
    for (t in 1:x$n_t){
      tryCatch(expr={
        if (x$pi[[i]][t]>0){
          cluster_it = clusteri %>% dplyr::filter(times==t)
          f1 <- approxfun(cluster_it$x, (cluster_it$y)/(binwidth*(x$pi[[i]][t])))  
          norm = integrate(f1, lower = 0.0, upper = 1, subdivisions=5000)$value
          f2 <- approxfun(cluster_it$x, (cluster_it$y)*(cluster_it$x)/(binwidth*(x$pi[[i]][t])))
          mean = integrate(f2, lower = 0.0, upper = 1, subdivisions=5000)$value
          f3 <- approxfun(cluster_it$x, (cluster_it$y)*(cluster_it$x - mean)**2/(binwidth*(x$pi[[i]][t])))
          var = integrate(f3, lower = 0.0, upper = 1, subdivisions=5000)$value
        } else {  mean = NA ; var = NA }
        mean_x[t] <- mean 
        var_x[t] <- var 
      }, error = function(e){
        mean_x[t] <- NA 
        var_x[t] <- NA
        }
      ) 
    }
    
    x$Clusters[which(x$Clusters$cluster == i & x$Clusters$type == "Mean"),]$fit.value <- list(mean_x)
    x$Clusters[which(x$Clusters$cluster == i & x$Clusters$type == "Variance"),]$fit.value <- list(var_x)
  }
  return(x)
}

get_permuted_order = function(x)
{ # permute beta clusters per timepoint in each possible way; find permutation that, per mutation, maximises product
  # of responsibilities over time; apply permutation to z_nkt and pdf.w
  
  Kbeta = x$Kbeta
  n_t = x$n_t
  
  # Create a sample 2D array
  matrix_data <- matrix(1:(Kbeta*n_t), nrow = Kbeta, ncol = n_t, byrow=TRUE)
  # Find all permutations within columns
  col_permutations <- apply(matrix_data, 2, FUN = function(col) permn(col))
  # Generate all combinations of col permutations
  combinations <- expand.grid(col_permutations)
  
  # separate from tail
  my_znkt = x$z_nkt[,2:x$K,]
  my_pdf.w = x$pdf.w[,2:x$K,]
  
  scores = rep(NA, length(combinations))
  
  for (i in 1:nrow(combinations)) {
    combination_matrix <- matrix(unlist(combinations[i, ]), nrow = nrow(matrix_data))
    
   # # TESTING -------------------------------
   # sorting_order = order(as.vector(t(combination_matrix)))
   # # to go from permuted to original
   # unsorted_combo = (as.vector(t(combination_matrix)))[sorting_order]
   # unsorted_combo = matrix(unsorted_combo,nrow = (K-1), ncol = n_t, byrow=TRUE)
   # # to go from original to permuted
   # sorted_orig = matrix(t(matrix_data)[as.vector(t(combination_matrix))],nrow = (K-1), ncol = n_t, byrow=TRUE)
   # #--------------------------------
      
    # permute z_nkt according to combination_matrix
    matlist_znkt = list()
    apply(my_znkt, MARGIN = c(1), FUN = function(x){
      y = matrix(t(x)[as.vector(t(combination_matrix))],ncol = Kbeta, nrow = n_t)
      #print(t(matrix(y[sorting_order],ncol = (K-1), nrow = n_t))) # to get back to original
      matlist_znkt[[length(matlist_znkt) + 1]] <<- t(y)
      })
    
    # convert back into 3D array
    znkt_permuted = matlist_znkt %>% abind(along=3) %>% aperm(c(3, 1, 2))
    
    # take product of permuted responsibilities across time
    resp_prods = apply(znkt_permuted, c(1,2), function(x){prod(x,na.rm = TRUE)})
    # sum across components and observations and save score
    scores[i] = sum(resp_prods,na.rm=TRUE)
  }
  my_permutation_index = which.max(scores)
  best_permutation = matrix(unlist(combinations[my_permutation_index, ]), nrow = nrow(matrix_data))
  
  # permute z_nkt and pdf.w according to best permutation
  matlist_znkt = list()
  matlist_pdf.w = list()
  apply(my_znkt, MARGIN = c(1), FUN = function(x){
    y = matrix(t(x)[as.vector(t(best_permutation))],ncol = Kbeta, nrow = n_t)
    matlist_znkt[[length(matlist_znkt) + 1]] <<- t(y)
    })
  apply(my_pdf.w, MARGIN = c(1), FUN = function(x){
    y = matrix(t(x)[as.vector(t(best_permutation))],ncol = Kbeta, nrow = n_t)
    matlist_pdf.w[[length(matlist_pdf.w) + 1]] <<- t(y)
    })
  znkt_permuted = matlist_znkt %>% abind(along=3) %>% aperm(c(3, 1, 2))
  pdf.w_permuted = matlist_pdf.w %>% abind(along=3) %>% aperm(c(3, 1, 2))
  
  # recombine with tail
  znkt_permuted = abind(x$z_nkt[,1,,drop=FALSE], znkt_permuted, along=2)
  pdf.w_permuted = abind(x$pdf.w[,1,,drop=FALSE], pdf.w_permuted, along=2)
  
  return(list(best.perm = best_permutation, z_nkt.perm = znkt_permuted, pdf.w.perm = pdf.w_permuted))
}  

# rename Beta clusters per timepoint according to permutation
rename_Beta_clusters_permutation = function(x)
{
  # get permutation
  perm_list = get_permuted_order(x) 
  my_permutation = perm_list$best.perm
  # replace x$z_nkt and x$pdf.w with permuted versions
  x$z_nkt = perm_list$z_nkt.perm
  x$pdf.w = perm_list$pdf.w.perm
  
  # permute Clusters object
  fit_values = x$Clusters %>% filter(cluster != 'Tail')
  tail_app = x$Clusters %>% filter(cluster == 'Tail')
  split_tib = split(fit_values, fit_values$type)
  permuted_clusters = as.list(rep(NA,length(names(split_tib))))
  i=0
  for (type in names(split_tib)){
   i=i+1
   split_tib[[type]] = split_tib[[type]] %>%
       separate(fit.value, into = paste0("fit.value_t_", 1:length(split_tib[[type]]$fit.value[[1]])), sep = ",") %>%
       separate(init.value, into = paste0("init.value_t_", 1:length(split_tib[[type]]$init.value[[1]])), sep = ",") %>%
       mutate(across((starts_with("fit.value_t_")|starts_with("init.value_t_")), ~gsub('[c\\(\\)\\" ]', '', .))) %>%
       mutate_all(~na_if(.x, "NA")) %>%
       mutate(across(starts_with("fit.value_t_")|starts_with("init.value_t_"), function(x){ifelse(!is.na(x),as.numeric(x),NA)}))
    
   a_df = split_tib[[type]] %>% select(-type)
   sub = a_df %>% arrange(cluster)
   
   fitval_sub = sub %>% select(-cluster) %>% select(starts_with("fit.value_t_")) 
   initval_sub = sub %>% select(-cluster) %>% select(starts_with("init.value_t_")) 
   
   fitval_permute = t(matrix(t(fitval_sub)[as.vector(t(my_permutation))],ncol = x$Kbeta, nrow = x$n_t)) %>% as_tibble()
   initval_permute = t(matrix(t(initval_sub)[as.vector(t(my_permutation))],ncol = x$Kbeta, nrow = x$n_t)) %>% as_tibble()
   
   names(fitval_permute) <- names(initval_permute) <- paste0('t_', seq(1, x$n_t)) 
   fitval_permute = fitval_permute %>% rowwise() %>% mutate(fit.value = list(c_across(starts_with('t_')))) %>% select(-starts_with('t_'))
   initval_permute = initval_permute %>% rowwise() %>% mutate(init.value = list(c_across(starts_with('t_')))) %>% select(-starts_with('t_'))
   permuted = list(fitval_permute,initval_permute) %>% bind_cols() %>% mutate(type = type) 
   permuted_clusters[[i]] = permuted %>% bind_cols(tibble(cluster = sub$cluster))
  } 
 tail_app = tail_app %>% mutate(init.value = lapply(init.value, function(x) list(x)))
 x$Clusters = bind_rows(permuted_clusters) %>% bind_rows(tail_app) %>% tibble::as_tibble()
 
 # setting x$a and x$b from new Clusters tibble 
 my_clus = x$Clusters %>% filter(type == 'a') %>% as_tibble() %>% select(-type,-init.value) 
 x$a =  my_clus %>% select(-cluster) %>% as.list(.$fit_value) %>% unname() %>% .[[1]]
 names(x$a) = my_clus$cluster
 my_clus = x$Clusters %>% filter(type == 'b') %>% as_tibble() %>% select(-type,-init.value) 
 x$b =  my_clus %>% select(-cluster) %>% as.list(.$fit_value) %>% unname() %>% .[[1]]
 names(x$b) = my_clus$cluster
 
 # this can be replaced by setting from Clusters as for a and b above
 # permute x$pi
 pi = x$pi %>% tibble::enframe(name='cluster',value='t') %>% unnest_wider(t, names_sep = "_")
 tail = pi %>% filter(cluster=='Tail')
 sub_pi = pi %>% arrange(cluster) %>% filter(cluster!='Tail') 
 sub_pi_permute = sub_pi %>% select(-cluster)
 col_names = names(sub_pi_permute)
 subpi_permute = t(matrix(t(sub_pi_permute)[as.vector(t(my_permutation))],ncol = x$Kbeta, nrow = x$n_t)) %>% as_tibble()
 names(subpi_permute) = col_names
 subpi_permute = subpi_permute %>% mutate(cluster = sub_pi$cluster)
 subpi_permute = bind_rows(tail, subpi_permute)
 x$pi = subpi_permute %>% rowwise() %>% 
   mutate(value = list(c_across(starts_with('t_')))) %>%
   select(-starts_with('t_')) %>%
   tibble::deframe()

 
 return(x)
}  

# Rename Beta clusters so that C1 is the one with highest mean etc.
rename_Beta_clusters = function(x)
{
  params = x$Clusters %>%
    filter(type == 'Mean', cluster != 'Tail') %>% select(-init.value) %>% 
    mutate(fit.value.mean = lapply(fit.value,mean))  %>%
    arrange(desc(fit.value.mean)) %>%
    mutate(new.name = paste0('C', row_number())) %>% select(-fit.value.mean)
  
  mapping = params$new.name
  names(mapping) = params$cluster
  
  mapping['Tail'] = 'Tail'
  mapping = mapping[c('Tail', params$cluster)]
  
  # Copy of x
  y = x
  
  # clusters
  y$data$cluster = mapping[y$data$cluster] %>% as.vector
  
  # numbers from clustering
  names(y$N.k) = mapping[names(y$N.k)] %>% as.vector
  
  # LV - forced assumed the order has been maintained
  colnames(y$z_nkt) =  mapping[colnames(y$z_nkt)]
  colnames(y$pdf.w) = mapping[colnames(y$pdf.w)]
  names(colnames(y$z_nkt)) = names(colnames(y$pdf.w)) = NULL 
  
  # Clusters table
  y$Clusters$cluster = mapping[y$Clusters$cluster]  %>% as.vector
  
  # mixing
  # names(y$pi) = mapping[names(y$pi)]  %>% as.vector
  y$pi = .params_Pi(y)
  
  # Beta parameters
  names(y$a) = mapping[names(y$a)]  %>% as.vector
  names(y$b) = mapping[names(y$b)]  %>% as.vector
  
  y
}

possible_NA_clusters <- function(n_t){
  # all combos 
  combos = c()
  if (n_t>1){
    for (i in 1:(n_t-1)){
      new_combos = combn(1:n_t, i)
      to_append = new_combos %>% asplit(MARGIN = 2)  %>% unname()
      combos = c(combos, to_append)
    } 
  }
  
  return(combos) 
}

possible_NA_clusters_4t <- function(t,n_t)
{ # finding all combinations of timepoints that include NA for separate post-hoc clustering
  combos = possible_NA_clusters(n_t) 
  
  # excluding those combos that include j
  if (length(combos)!=0){
    filtered_list <- lapply(combos, function(vec) { if (!(t %in% vec)) { return(vec) }})
    filtered_list <- filtered_list[sapply(filtered_list, function(x) !is.null(x))]
  } else {
    filtered_list = NULL
  }
  
  return(filtered_list)
}  

posthoc_NA_clusters <- function(data, offset, mask=FALSE){
  # offset is 0 if BL is first timepoint, 1 if next possible timepoint, etc.
  # mask: masking dataframe for use with hypothesis testing.
  #       Dataframe element interpretation:  
  #       TRUE  => null hypothesis of clonality can be rejected, mutation can 
  #                be considered truly missing
  #       FALSE => mutation can't be considered truly missing
  # Set mask parameter to TRUE for incorporation of hypothesis testing into 
  # calculation of mutation clusters
  n_t = mobster:::get_no_timepoints(data)
  
  data$cluster_orig = data$cluster
  data$clusterNA = data$cluster
  poss_NA_clusters = mobster:::possible_NA_clusters(n_t) 
  
  for (cl in poss_NA_clusters){
    cluster_name = paste('NA_in', paste(cl+offset, collapse = "_"), sep = "_")
    na_rows_bool = data$VAF %>% dplyr::select(as.character(cl)) %>%
      apply(., 1, function(row) all(is.na(row)))
    if ((mask==TRUE)&('sig_missing' %in% names(data))){
      pass_sig = data$sig_missing %>% 
        dplyr::select(paste0('X',as.character(cl))) %>%
        apply(., 1, function(row) all(row))
    } else {pass_sig = TRUE} # no hypothesis testing
    na_rows = which((na_rows_bool & pass_sig))
    
    data[na_rows,] = data[na_rows,] %>% dplyr::mutate(clusterNA = cluster_name)
  }
  # overwrite with NA if all timepoints are NA
  row_all_na <- apply(data$VAF, 1, function(row) all(is.na(row)))
  data[row_all_na,] = data[row_all_na,] %>% mutate(clusterNA = NA)
  
  return(data)
}

# create tibble for assignment of initial kmeans clustering 
.cluster_kmeans = function(fit){
  kmeans_clusters = tibble()
  for (t in 1:fit$n_t){
    new_col = fit$z_nkt[,,t] %>% apply(MARGIN=1,
               function(x) {names(fit$pi)[which(x == max(x, na.rm = TRUE))[1]]})
    if (t==1){
      kmeans_clusters = as_tibble(new_col)
    } else {
      kmeans_clusters = bind_cols(kmeans_clusters, new_col)
    }
  }
  names(kmeans_clusters) = names(fit$data$VAF)
  return(kmeans_clusters)
}

