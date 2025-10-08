
time_in_units <- function(observation_days){
  # second observation is at t=0 (for increased symmetry)
  nos_days = lubridate::interval(observation_days[2],observation_days)/ddays(1)
  my_times = nos_days/30 # 1 unit = 30 days
  return(my_times)
}

get_no_timepoints <- function(df){return(length(names(df$VAF)))}

calc_poly_t <- function(coeffs,t){
  val = 0
  for (i in 1:length(coeffs)){
    val = val + coeffs[i]*t**(i-1)
  }
  return(val)
}

calc_poly_t_all <- function(coeffs,t){
  vals = sapply(t,function(x){calc_poly_t(unlist(coeffs),x)})
  return(vals)
}

