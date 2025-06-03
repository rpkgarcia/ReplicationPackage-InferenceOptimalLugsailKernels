
AR1_data <- function(rho_y=.7,  e_sd=1, big_T=200, mean_vec = 0){
  
  # Generate the data.  Initial value is 0
  sim_data <- rep(0, big_T)
  
  # The rest of the values
  for(t in 2:big_T){
    e_t <- rnorm(1, 0, e_sd)
    
    sim_data[t] <- mean_vec + rho_y*sim_data[c(t-1)] + e_t
  }
  
  return(sim_data)
}
