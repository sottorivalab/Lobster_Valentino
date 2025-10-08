
# expected VAF of H0 prob. distribution according to multiplicity m
VAF_m <- function(m, purity, ccf, ploidy){
  vaf_m = (m*purity*ccf)/(2*(1-purity) + purity*ploidy) # from https://doi.org/10.1186/s13059-024-03170-5
  return(vaf_m)
}

# gives prob. mass of binomial distribution at observed count 0 for multiplicity m
elementwise_binom_prob0 <- function(depth, m, purity, ccf, ploidy) { # m = multiplicity
  return(dbinom(x=0, size = depth, prob = min(1,VAF_m(m, purity, ccf, ploidy)) ))  
}

elementwise_sig_missing <- function(x, y, correction=1) {
   return(ifelse(is.na(x),y<(0.05/correction), NA))  # For example, you can use any function here
}

calc_binomial_weights <- function(ploidy){ 
  
  if (is.na(ploidy)){ ploidy=2 } # assume NA follows diploid
  weight_nums = c()
  multiplicities = 1:(ploidy-1)
  
  if (ploidy!=1){
    sum = 0
    for (i in multiplicities){
      weight_num = dbinom(i,ploidy,0.5)
      sum = sum + weight_num
      weight_nums = c(weight_nums,weight_num)
    }   
  } else { # monoploid is special case
    weight_num = 1
    sum = 1
    weight_nums = weight_num
  }
  multiplicities = multiplicities[multiplicities!=0]
  
  for (i in multiplicities){
    weight_nums[i] = weight_nums[i]/sum
  }
  
  return(list(m = multiplicities, w = weight_nums))
}
  
# weights for binomial mixture of binomial distributions for different ploidies, excepting 0 and 1
binomial_mixture_probdist = function(ploidy){
  
  if (is.na(ploidy)){ ploidy=2 } # assume NA follows diploid
  
  func_list = list()
  
  cbw = calc_binomial_weights(ploidy)
  multiplicities = cbw$m
  weight_nums = cbw$w
  
  for (m in 1:length(multiplicities)){
    if (m!=0){
      func_list[[m]] = substitute(function(depth, purity, ccf, ploidy){
         return( weight * elementwise_binom_prob0(depth,const,purity,ccf,ploidy) ) }, list(weight = weight_nums[m], const = m) )
    }
  } # apply eval() to get function
  
  return(func_list)
} 

# evaluate funclist for a certain timepoint
eval_funclist_elementwise = function(func_list, depth, purity, CCF_H0, ploidy){
  result = 0
  # just a sum over all elements of func_list
  for (i in 1:length(func_list)){
     result = result + eval(func_list[[i]])(depth, purity, CCF_H0, ploidy)  
  }
  return(result)
}

# Wrapper function to apply element-wise operation and handle NA values
safe_eval_funclist_elementwise <- function(func_list, depth, purity, CCF_H0, ploidy) {
  if (is.na(depth) || is.na(purity) || is.na(CCF_H0) || is.na(ploidy) || is.null(func_list)) {
    return(NA)
  } else {
    return(eval_funclist_elementwise(func_list, depth, purity, CCF_H0, ploidy))
  }
}
     
## calculating percentiles using purities as weights (excluding NAs)
#calc_percentile <- function(df, purities){
#  df$percentile = df$VAF_orig
#  df$percentile <- apply(df$percentile, 2, function(x) {rank(x,na.last = "keep") / (sum(!is.na(x)) + 1)}) %>%
#    as.data.frame()
#  df$weighted_percentile = sapply(1:nrow(df$percentile), function(x) {weighted.mean(df$percentile[x,], purities, na.rm = TRUE)})
#  
#  return(df)
#}

# calculating percentiles using sequencing depths as weights (excluding NAs)
calc_percentile <- function(df){
  df$percentile = df$VAF_orig
  df$percentile <- apply(df$percentile, 2, function(x) {rank(x,na.last = "keep") / (sum(!is.na(x)))}) %>%
    as.data.frame()
  df$weighted_percentile = sapply(1:nrow(df$percentile), function(x) {weighted.mean(df$percentile[x,], df$DP[x,], na.rm = TRUE)})
  
  return(df)
}

# Benjamini-Hochberg procedure
# find significance (whether null hypothesis can be rejected at level alpha) using Benjamini-Hochberg procedure
sig_BH_proc <- function(VAF_orig, prob_H0, alpha=0.05){
  
  flatten_df_as_vec <- function(df){
    df_asvec = unlist(df)
    names(df_asvec) = df_asvec
    return(df_asvec)
  }
  
  # flatten dataframes and order in terms of ascending prob_H0 
  serialised_prob_H0 = data.frame(Values = flatten_df_as_vec(prob_H0))   
  permutation = order(flatten_df_as_vec(serialised_prob_H0))
  serialised_prob_H0 = serialised_prob_H0[permutation,]
  serialised_VAF_orig = data.frame(Values = flatten_df_as_vec(VAF_orig))
  serialised_VAF_orig = serialised_VAF_orig[permutation,]
  
  reduced_prob_H0 = serialised_prob_H0[is.na(serialised_VAF_orig)]
  m = length(reduced_prob_H0)
  # finding largest k such that prob_H0(k) <= k*alpha/m
  thresholds = ((1:m)/m)*alpha
  k=1
  H0_rejections = c()
  for (i in 1:length(serialised_prob_H0)){
    H0_rejection = NA
    if(is.na(serialised_VAF_orig[i])){
      H0_rejection = serialised_prob_H0[i]<=thresholds[k]
      k = k+1
    }
    H0_rejections = c(H0_rejections,H0_rejection)
  }
  
  # reversing permutation
  reverse_perm = match(1:length(permutation), permutation) 
  H0_rejections = H0_rejections[reverse_perm] 
    
  # reforming dataframe
  parts <- split(H0_rejections, ceiling( seq_along(H0_rejections) / ceiling(length(H0_rejections)/ncol(prob_H0)) ))
  H0_rejections = data.frame(parts)
  
  return(H0_rejections)
}

# calculate VAF at highest peak
calc_VAF_clone <- function(fit){
  if (is.null(fit$best$hc_attribs)){ fit$best = hard_clustered_VAF_CCF_tail(fit$best) }
  means = fit$best$hc_attribs %>% dplyr::filter(type=='VAF_mean') %>% select(value)
  vectors <- means[[1]]
  vectors = lapply(vectors,function(x){replace(x, is.na(x), -Inf)}) # replace na values with -Inf
  if (length(vectors)>1){
    VAF_clone <- do.call(pmax, vectors) # for now, just ignore NA
  } else { VAF_clone = unlist(vectors) }
  
  return(VAF_clone)
}
          
# calculating expected VAFs from weighted percentile and observed VAF distribution for each ploidy fit
calc_CCF_H0 <- function(ploidy_fit){
  
  calc_VAF_expected_t <- function(t, fit){
    
    w_percentile = fit$weighted_percentile
    VAF = fit$VAF_orig[t]
    VAF_sans_NA = unlist(VAF)[which(!is.na(VAF))]
    
    expected_VAF = rep(NA,length(w_percentile))
    
    for (i in 1:length(w_percentile)){
      
      if (!is.na(fit$VAF_orig[i,t])){ 
        # mutation not missing
        expected_VAF[i] = NA 
      } else {
        value <- w_percentile[i]
        
        # find quantile of percentile
        quantiles <- quantile(VAF_sans_NA, probs = value)
        # closest element 
        if (!is.na(quantiles)|!is.na(quantiles)){
          expected_VAF[i] <- VAF_sans_NA[which.min(abs(VAF_sans_NA - quantiles))] 
        } else {
          expected_VAF[i] = NA 
        }
      }
    }
    
    return(expected_VAF)
  }
  
  calc_CCF_expected_t <- function(t, fit){
    VAF_expected_t = calc_VAF_expected_t(t, fit$annotated) 
    VAF_clone_t = calc_VAF_clone(fit)[t]
    CCF_expected_t = VAF_expected_t/VAF_clone_t
    
    return(CCF_expected_t)
  }
  
  CCF_H0 <- sapply(1:ncol(ploidy_fit$annotated$VAF_orig), function(x){calc_CCF_expected_t(x,ploidy_fit)}) %>% as.data.frame()
  
  return(CCF_H0) 
}