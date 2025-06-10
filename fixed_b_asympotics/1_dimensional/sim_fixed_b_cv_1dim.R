# -------------------------------------------------------------------------
# Load Functions ----------------------------------------------------------
# -------------------------------------------------------------------------


# Bartlett 
bartlett <- function(x){
  if(abs(x)<1){
    k_x <- 1-abs(x)
  } else{
    k_x <- 0 
  }
  return(k_x)
}

# Quadratic Spectral 
qs <- function(x){
  p1 <- sin(6*pi*x/5)/(6*pi*x/5)
  p2 <- cos(6*pi*x/5)
  p3 <- 25/(12*pi^2*x^2)
  k_x <- p3*(p1-p2) 
  if(x == 0){
    k_x <- 1
  }
  return(k_x)
}

# Parzen
parzen <- function(x){
  k_x <- 0 
  if(abs(x)<=1){
    if(abs(x)<= 0.5){
      k_x <- 1 - 6*x^2 + 6*abs(x)^3
    } else {
      k_x <- 2*(1-abs(x))^3
    }
  } 
  return(k_x)
}



# Lugsail Transformation (main)
lugsail <- function(x, lugsail_parameters, the_kernel= bartlett){
  r <- lugsail_parameters$r 
  c <- lugsail_parameters$c
  
  # Actual lugsail 
  y1 <- the_kernel(x)/(1-c) 
  y2 <- 0 
  
  # If QS, then always include the extra bit. 
  if(deparse(substitute(the_kernel)) == "qs"){
    y2 <- the_kernel(x*r)*c/(1-c) 
  }
  
  # If not QS, then you need to check
  if(abs(x) < 1/r){
    y2 <- the_kernel(x*r)*c/(1-c) 
  }
  y <- y1- y2
  
  return(y)
}

# Lugsail Support Function to get lugsail parmeters  
# (default) b = Andrews (1991) Rule: 0.75*big_T^(-2*q/(2*q+1))
get_lugsail_parameters <- function(big_T, q, method = "Zero", 
                                   b = 0.75*big_T^(-2*q/(2*q+1))){
  
  if(method == "Over"){
    r <- 3
    c <- 2/(1+r^q)
    
  } else if(method == "Adaptive"){
    r <- 2
    M  <- big_T * b
    c_num <- (log(big_T) - log(M) + 1)
    c_den <- r^q*(log(big_T) - log(M)) + 1
    c <- c_num/c_den 
    
  } else {
    # Zero or Manual lugsail
    r <- 2
    c <- r^(-q)
    
  }
  parameters <- list(r = r, c = round(c, 2))
  return(parameters)
}

# Calculate all the autocovariances for a 1-dimensional simulation
all_R <- function(one_sim_data){
  big_T <- length(one_sim_data)
  all_auto_cov <- sapply(0:(big_T-1), R, one_sim_data= one_sim_data)
  return(all_auto_cov)
}



# one_sim_data = a vector with one simulation
R <- function(h, one_sim_data){
  big_T <- length(one_sim_data)
  est_mean <- mean(one_sim_data)
  index <- 1:(big_T -h)
  est <- (one_sim_data[index]-est_mean)*(one_sim_data[(index+h)] - est_mean)
  
  autocov_s <- sum(est)/big_T
  
  if(h!=0){
    autocov_s <- 2*autocov_s
  }
  
  return(autocov_s)
}


# -------------------------------------------------------------------------
# Estimate LRV  --------------------------------------------------
# -------------------------------------------------------------------------

LRV_estimator <- function(try_b, all_sim_data, the_means, all_autocovariances, 
                          the_kernel, lugsail_parameters = list(r = 1, c= 0), 
                          mother_omega){
  # Make the wieghts that correspond to the autocovariances
  # Each row corresponds to the weights for a specific weight. 
  # Each column corresponds to the autocovariance lag
  # Each simulation gets each row of these scenarios. 
  W <- sapply(try_b, function(b){
    M <- b*big_T 
    new_weights <- sapply(0:(big_T -1)/(M), lugsail, the_kernel = the_kernel, 
                          lugsail_parameters = lugsail_parameters)
  })
  W <- t(W)
  rownames(W) = paste("b=", round(try_b, 3), sep = "")
  colnames(W) = paste("Lag=", 0:(nrow(all_sim_data)-1), sep = "")
  
  
  omega_hats <- W %*% all_autocovariances
  colnames(omega_hats) <- paste("Sim", 1:ncol(all_sim_data))
  rownames(omega_hats) <- paste("b=", round(try_b, 3))
  
  # Fix PSD Issues 
  index <- (omega_hats <= 0)
  omega_hats[index] <- mother_omega[index]
  
  # Return Statistics 
  psd_counts <- apply(index, 1, sum)
  
  return(list(omega_hats, psd_counts))
}



LRV_mother_estimator <- function(try_b, all_sim_data, the_means,
                                 all_autocovariances, the_kernel){
  # Make the wieghts that correspond to the autocovariances
  # Each row corresponds to the wieghts for a specific wieght. 
  # Each column corresponds to the autocovariance lag
  # Each simulation gets each row of these scenarios. 
  W <- sapply(try_b, function(b){
    M <- b*big_T 
    new_weights <- sapply(0:(big_T -1)/(M), the_kernel)
  })
  W <- t(W)
  rownames(W) = paste("b=", round(try_b, 3), sep = "")
  colnames(W) = paste("Lag=", 0:(nrow(all_sim_data)-1), sep = "")
  
  
  omega_hats <- W %*% all_autocovariances
  colnames(omega_hats) <- paste("Sim", 1:ncol(all_sim_data))
  rownames(omega_hats) <- paste("b=", round(try_b, 3))
  return(omega_hats)
}

# -------------------------------------------------------------------------
# Calculate F-Statistics --------------------------------------------------
# -------------------------------------------------------------------------

# Calculate the F test statistic for dimensions 1, ..., d_max
# Under the 

# b: the bandwidth, proportion of estimated autocovs that have non-zero weight
# the_sim_data: 
#     - Contains `big_T` simulated dependent random vectors of dimension (d_max).
#     - Should already be centered with hypothesized or estimated means.
# the_means: 
#     - The estimated mean vector.  
#     - Should be of length d.
#     - Need this because the_sim_data is already centered. 
# all_autocovariances: 
#     - Rows are correspond to estimated autocov (R) at lag [0, ..., (big_T -1)]
#     - Columns correspond to the vectorization of the estimated autocov matrix: 
#          R11, R12, R13, ..., R1d, R21, R22, ..., R2d, ..., Rd1, ...Rdd
# kernel: the name of the kernel function 
# lugsail_parameters: 
#     - a named list, list(r = 1, c= 0) that contains the lugsail parameters
#     - default is non-lugsail  

F_stats <- function(the_means, omega_hats, null_mean = 0){
  
  # ------- F-statistics for various b values -------
  # [ #,  ] a different b value
  # [  , #] a different simulation
  omega_hat_inv <- 1/omega_hats
  num <- (the_means- null_mean)
  F_stat <- apply(omega_hat_inv, 1, function(row) row*num^2*big_T)
  F_stat <- t(F_stat)
  
  return(F_stat)
}


# -------------------------------------------------------------------------
# Main --------------------------------------------------------------------
# -------------------------------------------------------------------------

get_kernel_F_stats <- function(try_b, all_sim_data, the_means,
                               all_autocovariances, the_kernel, q = 1){
  
  # Mother Kernel 
  omega_mother <- LRV_mother_estimator(try_b, all_sim_data, the_means,
                                       all_autocovariances, the_kernel)
  F_mother <- t(F_stats(the_means, omega_mother))
  
  # Zero Lugsail 
  lug_para <- get_lugsail_parameters(big_T, q = q, method = "Zero")
  omega_zero <- LRV_estimator(try_b, all_sim_data, 
                              the_means, all_autocovariances, 
                              the_kernel = the_kernel, 
                              lugsail_parameters = lug_para, 
                              omega_mother)
  F_zero <- t(F_stats(the_means, omega_zero[[1]]))
  
  # Adaptive Lugsail 
  big_T <- nrow(all_sim_data)
  init_b <-  0.75*big_T^(-2*q/(2*q+1))
  lug_para <- get_lugsail_parameters(big_T, q = q, method = "Adaptive", 
                                     b = init_b)
  omega_adapt <- LRV_estimator(try_b, all_sim_data, 
                               the_means, all_autocovariances, 
                               the_kernel= the_kernel, 
                               lugsail_parameters = lug_para, 
                               omega_mother)
  F_adapt <- t(F_stats(the_means, omega_adapt[[1]]))
  
  # Over Lugsail 
  lug_para <- get_lugsail_parameters(big_T, q = q, method = "Over")
  omega_over <- LRV_estimator(try_b, all_sim_data, 
                              the_means, all_autocovariances, 
                              the_kernel= the_kernel,
                              lugsail_parameters = lug_para, 
                              omega_mother)
  F_over <- t(F_stats(the_means, omega_over[[1]]))
  
  return(list(F_stats = list(mother = F_mother, zero = F_zero, 
                             adapt = F_adapt, over = F_over), 
              psd_counts = data.frame(zero = omega_zero[[2]], 
                                      adapt = omega_adapt[[2]], 
                                      over = omega_over[[2]])))
}


# -------------------------------------------------------------------------
# Main --------------------------------------------------------------------
# -------------------------------------------------------------------------


# Generates a null data set.  Calculates the autocovariance matrices. 
# Calculates the F-statistic for dimensions (1,..., d_max). 
# Need to repeat this simulation many times to then find what asymptotic CV is. 


simulate_f_stat <- function(big_T = 1000, num_sims = 100, null_mean = 0){
  
  # ------- Simulate all of the data  -------
  # [#,  ] each value for a single simulation, 1:big_T
  # [ , #]  each simulation 1:num_sims 
  all_sim_data <- replicate(num_sims, rnorm(big_T))
  orig_sim_data <- all_sim_data
  
  # Location model, get error terms 
  the_means <- colMeans(all_sim_data)
  all_sim_data <- all_sim_data- the_means
  
  # ------- AutoCovariance Matrices  -------
  # [#,  ] the lag, (0, ..., big_T-1)
  # [ , #]  each simulation, 1:num_sims 
  all_autocovariances <- apply(all_sim_data, 2, all_R)
  
  # ------- F-statistics for various b values -------
  # [ #,  ] a different b value
  # [  , #] a different simulation
  
  # ------- BARTLETT ------- 
  F_bartlett <- get_kernel_F_stats(try_b, all_sim_data, the_means,
                                   all_autocovariances, bartlett, q = 1)
  psd_bartlett <- F_bartlett[[2]]
  F_bartlett <- F_bartlett[[1]]
  
  #  ------- PARZEN -------
  F_parzen <- get_kernel_F_stats(try_b, all_sim_data, the_means,
                                 all_autocovariances, parzen, q = 2)
  psd_parzen <- F_parzen[[2]]
  F_parzen <- F_parzen[[1]]
  
  # -------  QS -------
  F_qs <- get_kernel_F_stats(try_b, all_sim_data, the_means,
                             all_autocovariances, qs, q = 2)
  psd_qs <- F_qs[[2]]
  F_qs <- F_qs[[1]]
  
  
  # ------- Store all F Stats -------
  # [ #,  ,  ] the simulation number 
  # [  , #,  ] a different b value
  # [  ,  , #] a the_kernel 
  F_stats_all <- array(c(F_bartlett[["mother"]], F_bartlett[["zero"]], 
                         F_bartlett[["adapt"]], F_bartlett[["over"]], 
                         F_parzen[["mother"]], F_parzen[["zero"]], 
                         F_parzen[["adapt"]], F_parzen[["over"]], 
                         F_qs[["mother"]], F_qs[["zero"]], 
                         F_qs[["adapt"]], F_qs[["over"]]),
                       dim=c(num_sims, length(try_b),  12))
  
  dimnames(F_stats_all)[[3]] <-  c("Bartlett_Mother", "Bartlett_Zero",
                                   "Bartlett_Adapt", "Bartlett_Over", 
                                   "Parzen_Mother", "Parzen_Zero",
                                   "Parzen_Adapt", "Parzen_Over", 
                                   "QS_Mother", "QS_Zero",
                                   "QS_Adapt", "QS_Over")
  
  
  # ------- Store PSD Statistics -------
  psd_counts <- list(bartlett = psd_bartlett, 
                     parzen = psd_parzen, 
                     qs = psd_qs)
  
  return(list(F_stats_all, psd_counts))
}

# -------------------------------------------------------------------------
# Simulation Parameters (Adjust me) ---------------------------------------
# -------------------------------------------------------------------------

#  What bs to use 
try_b <-  seq(0.005, .99, by = 0.005)

# How many replicates
# KV005 used 50,000
num_sims <- 50000

# Sample size of each replicate
big_T = 1000


# -------------------------------------------------------------------------
# Run Simulation  ---------------------------------------------------------
# -------------------------------------------------------------------------
set.seed(62)
all_F_stats <- simulate_f_stat(big_T = 1000, num_sims = num_sims, null_mean = 0)
psd_counts <- all_F_stats[[2]]
all_F_stats <- all_F_stats[[1]]

# [ #,  ,  ] : each simulation [1:num_sims]
# [  , #,  ] :rows, the different b values [1:length(try_b)]
# [  ,  , #] : the original kernel (1), and lugsail (2) [1:2]

dimnames(all_F_stats)[[1]] <- paste("sim", 1:num_sims, sep="")
dimnames(all_F_stats)[[2]] <- round(try_b,3)
dimnames(all_F_stats)[[3]] <- c("Bartlett_Mother", "Bartlett_Zero",
                                "Bartlett_Adapt", "Bartlett_Over", 
                                "Parzen_Mother", "Parzen_Zero",
                                "Parzen_Adapt", "Parzen_Over", 
                                "QS_Mother", "QS_Zero",
                                "QS_Adapt", "QS_Over")

# -------------------------------------------------------------------------
# Create CV Tables --------------------------------------------------------
# -------------------------------------------------------------------------

save_cv <- function(the_kernel){
  F_stats <- all_F_stats[ , ,the_kernel]
  
  t90 <- apply(F_stats, 2, quantile, probs = 0.90)
  write.csv(t90, paste(the_kernel, "_10.csv", sep = ""))
  
  t95 <- apply(F_stats, 2, quantile, probs = 0.95)
  write.csv(t95, paste(the_kernel, "_05.csv", sep = ""))
  
  t97.5 <- apply(F_stats, 2, quantile, probs = 0.975)
  write.csv(t97.5, paste(the_kernel, "_025.csv", sep = ""))
  
  t99 <- apply(F_stats, 2, quantile, probs = 0.99)
  write.csv(t99, paste(the_kernel, "_01.csv", sep = ""))
}

#  this is how I get all the critical values 
setwd("Fixed_b_CV_Tables")
sapply(dimnames(all_F_stats)[[3]], save_cv)
setwd("..")


# -------------------------------------------------------------------------
# Save for Distribution Plots ---------------------------------------------
# -------------------------------------------------------------------------

# Saves a data frame for a specific the_kernel and dimension. 
# Keeps the F-statistics so I can use them to make distribution plots later. 

save_distr <- function(the_kernel){
  F_stats <- all_F_stats[ , ,the_kernel]
  write.csv(F_stats, paste(the_kernel,".csv", sep = ""))
}

#  this is how I get all the critical values 
setwd("Fixed_b_distribution")
sapply(dimnames(all_F_stats)[[3]], save_distr)
setwd("..")



# -------------------------------------------------------------------------
# PSD Plots ---------------------------------------------
# -------------------------------------------------------------------------

# Saves a data frame for a specific the_kernel and dimension. 
# Keeps the F-statistics so I can use them to make distribution plots later. 

save_psd_info <- function(the_kernel){
  psd_info <- psd_counts[[the_kernel]]/num_sims
  write.csv(psd_info, paste(the_kernel,".csv", sep = ""))
}

#  this is how I get psd rate_plots 
setwd("psd_issues")
sapply(names(psd_counts), save_psd_info)
setwd("..")

