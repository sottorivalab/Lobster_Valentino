get_a_b <- function(data,des_min=0.001,des_max=0.99){
  # need to input ln of all VAFs as data variable
  dmin = min(data,na.rm=TRUE)
  dmax = max(data,na.rm=TRUE)
  a = (des_max-des_min)/(dmax-dmin)
  b = des_min - dmin*(des_max-des_min)/(dmax-dmin)
  return(list(a = a, b = b))
}

.custom_MinMaxScale <- function(data,des_min=0.001,des_max=0.99){
  # des_min/max: desired range values 
  dmin = min(data,na.rm=TRUE)
  dmax = max(data,na.rm=TRUE)
  scaled = ((data-dmin)/(dmax-dmin))*(des_max-des_min) + des_min
  return(scaled)
}

.cMMS_inverter <- function(data,des_min,des_max,dmin,dmax){
  # des_min/max: desired range values 
  # dmin/dmax: min/max values of original unscaled VAF
  scaled = ((data-des_min)/(des_max-des_min))*(dmax-dmin) + dmin
  return(scaled)
}

.log_untransform <- function(val,df){
  # undo x log transform and min-max scaling
  # error here
  un_tfd = exp(mobster:::.cMMS_inverter(val,
                                        min(df$VAF_ln_scaled,na.rm=TRUE),
                                        max(df$VAF_ln_scaled,na.rm=TRUE),
                                        min(df$VAF_ln,na.rm=TRUE),
                                        max(df$VAF_ln,na.rm=TRUE)))
  return(un_tfd)
}

.log_transform <- function(df){
  # do x log transform and min-max scaling
  new_cols = df$VAF %>% mutate_all(~log(.))
  colnames(new_cols) <- as.character(1:get_no_timepoints(df))
  MM_scaled = .custom_MinMaxScale(new_cols) %>% tibble()
  df = df %>% mutate("VAF_ln"= new_cols) %>% mutate("VAF_ln_scaled"= MM_scaled) %>% 
    mutate(VAF_orig = VAF) %>% mutate(VAF = VAF_ln_scaled) 
  return(df)
}


.log_transform_a_b <- function(df,a,b){
  new_cols = df$VAF %>% mutate_all(~log(.))
  MM_scaled = new_cols %>% mutate_all(~(a*(.) + b))
  colnames(new_cols) <- as.character(1:get_no_timepoints(df))
  df = df %>% mutate("VAF_ln"= new_cols) %>% mutate("VAF_ln_scaled"= MM_scaled) %>% 
    mutate(VAF_orig = VAF) %>% mutate(VAF = VAF_ln_scaled) 
  return(df)
}
