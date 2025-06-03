#setwd("~/Documents/GitHub/DissertationSim/simulations/CompareProcedures")
library(Matrix)
library(MASS)
library(distr)
D <- DExp(rate = 1) 



# Model Simulation --------------------------------------------------------

# Call outside of function for computation time. 
# exponent <- abs(matrix(1:10 - 1, nrow = 10, ncol = 10, byrow = TRUE)-(1:10 - 1))
# Sigma <- 0.9^exponent

AR_X <- function(big_T, rho, d){
  X <- matrix(0, nrow = big_T, ncol = d)
  
  for(i in 2:big_T){
    X[i, ] <- rho*X[(i-1), ] + mvrnorm(n = 1, rep(0, d), diag(d)) 
  }
  return(X)
}

AR_u <- function(big_T, rho){
  u <- rep(0, big_T)
  for(i in 2:big_T){
    u[i] <- rho*u[(i-1)] + rnorm(1)
  }
  return(u)
}


AR1_HOMO <- function(big_T, rho, d, theta = rep(0, d)){
  X <- AR_X(big_T, rho, d)
  u <- AR_u(big_T, rho)
  Y <- X%*%theta + u
  return(list(Y = c(Y), X = X))
}



HET_u <- function(big_T, rho, a_0 = 5, a_1 = .25){
  u <- rep(0, big_T)
  v <- AR_u(big_T, rho)
  
  for(i in 2:big_T){
    u[i] <- sqrt(a_0 + a_1*u[(i-1)])*v[i]
  }
  return(u)
}

AR1_HET <- function(big_T, rho, d, theta = rep(0, 4)){
  X <- AR_X(big_T, rho, d)
  u <- HET_u(big_T, rho)
  Y <- X%*%theta + u
  return(list(Y = c(Y), X = X))
}


ARMA_Gu <- function(big_T, rho, a = 0.5){
  Z <- rnorm(big_T)
  u <- rep(0, big_T)
  
  for(i in 2:big_T){
    u[i] <- rho*u[(i-1)] + a*Z[(i-1)] + Z[i] + rnorm(1)
  }
  return(u)
}

ARMA_G <- function(big_T, rho, d, theta = rep(0, 4)){
  X <- AR_X(big_T, rho, d)
  u <- ARMA_Gu(big_T, rho)
  Y <- X%*%theta + u
  return(list(Y = c(Y), X = X))
}



ARMA_Lu <- function(big_T, rho, a = 0.5){
  Z <- rnorm(big_T)
  u <- rep(0, big_T)
  
  for(i in 2:big_T){
    u[i] <- rho*u[(i-1)] + a*Z[(i-1)] + Z[i] + r(D)(1)
  }
  return(u)
}

# ARMA(1,1) with laplace disturbance 
ARMA_L <- function(big_T, rho, d, theta = rep(0, 4)){
  X <- AR_X(big_T, rho, d)
  u <- ARMA_Lu(big_T, rho)
  Y <- X%*%theta + u
  return(list(Y = c(Y), X = X))
}


# Practice calling 
errors_and_thetas <- function(big_T, rho, model_type, d, theta){
  if(model_type == "AR1_HOMO"){
    the_data <- AR1_HOMO(big_T, rho, d, theta)
  } else if(model_type == "AR1_HET"){
    the_data <- AR1_HET(big_T, rho, d, theta)
  } else if(model_type == "ARMA_G"){
    the_data <- ARMA_G(big_T, rho, d, theta)
  } else if(model_type == "ARMA_L"){
    the_data <- ARMA_L(big_T, rho, d, theta)
  } 
  
  X <- cbind(1, the_data$X)
  # fit <- lm(the_data$Y ~ the_data$X) # To check work (good)
  #X <- rep(1, big_T)
  coefs <- (solve(t(X)%*%X)%*% t(X)%*%the_data$Y)  # All coefs 
  M <- t(X)%*%X/big_T    # Do not consider intercept
  errors <- apply(X, 2, function(col) col*c(the_data$Y - X%*%coefs))  # errors to make Omega, not normal errors
  errors <- errors - apply(errors, 2, mean)  # Center them immediately 
  
  return(list(coefs= coefs[2], errors = errors, M = M))
}

# CV_Grid_Look_Up ---------------------------------------------------------


grid_cv_find <- function(d, alpha, the_kernel, lug_type, new_b){
  my_dir <- getwd()
  # Go to correct dimension 
  # setwd("../../fixed_b_asymptotics/1_dimensional/Fixed_b_CV_tables_New")
  
  # Make alpha a character
  alpha <- alpha*100
  if(alpha == 10){
    alpha <- "10"
  } else {
    alpha <- paste("0", alpha, sep = "")
  }
  
  file <- paste(the_kernel, lug_type, alpha, sep = "_")
  grid <- read.csv(paste("../../fixed_b_asymptotics/1_dimensional/Fixed_b_CV_Tables_New/", file, ".csv", sep = ""))
  
  cv_by_b <- sapply(new_b, function(b){
    return(grid[which.min(abs(b-grid[,1])),2])
  })
  #setwd(my_dir)
  
  return(cv_by_b)
}


# b Rules -----------------------------------------------------------------
mse_rule <- function(big_T, rho = .25, q = 1, d = 1){
  #the_b <- 0.75*big_T^(-2*q/(2*q+1))
  w_1 <- (2*rho/(1-rho^2))
  the_b  <- 1.1447*(w_1/big_T)^(2/3)
  if(is.na(the_b)){the_b <- 0}
  return(the_b)
}


# Flat Top Rule from Politis 2003
ft_rule <- function(all_correlation, big_T, c = 2){
  K_T <- max(5, floor(log(big_T)))
  bound <- c*sqrt(log(big_T)/big_T)
  k <- 1:K_T
  m_max <- length(all_correlation)-1
  for(m in 1:m_max){
    indices <- m+k 
    indices <- indices[indices<m_max]
    if(all(abs(all_correlation[indices])<bound)){
      break 
    }
  }
  m <- max(1, m-1)  
  the_b <- (2*m)/big_T
  return(the_b)
}

# Testing Optimal Rule (Mine) 
opt_rule <- function(rho, big_T, alpha = 0.05, d =1, lug_type = "Zero"){
  cv <- qchisq((1-alpha), d)/d
  num2 <- 0 
  w_q <- 2*rho/(1-rho^2)
  
  if (lug_type == "Over"){
    b_init <- opt_rule(rho, big_T, alpha, d, "Zero")
    g_q <- -1
    num2 <-  dchisq(cv, d)*cv*g_q*w_q/(b_init*big_T)
  }
  
  num1 <-  alpha/(big_T*log(rho))
  num1 <-  alpha^.5/(big_T*log(rho))
  den1 <- dchisq(cv, d)*cv*2*rho^2/(1+ rho)
  the_b <- log(-(num1 + num2)/den1)/(big_T*log(rho))
  the_b <- max(0, the_b)
  if(is.na(the_b)){the_b <- 0 }
  
  return(the_b)
}

# opt_rule(rho = .7, big_T = 200, alpha = 0.05, d = 1, lug_type = "Mother")
# opt_rule(rho = .7, big_T = 200, alpha = 0.05, d = 1, lug_type = "Zero")
# opt_rule(rho = .7, big_T = 200, alpha = 0.05, d = 1, lug_type = "Over")

# Kernel Functions --------------------------------------------------------



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



# source(paste(dir, "est_autocov.R", sep = ""))
R <- function(h, the_sim_data){
  big_T <- nrow(the_sim_data)
  index <- 1:(big_T -h)
  
  # Already centered
  est <- lapply(index, function(i){
    est <- the_sim_data[i, ]%*%t(the_sim_data[(i+h), ])
  })
  
  # Sum together
  autocov_s <- Reduce('+', est)/big_T
  
  # Because of symmetry 
  if(h!=0){
    autocov_s <- autocov_s + t(autocov_s)
  }
  
  return(autocov_s)
}


# LRV Estimator -----------------------------------------------------------


# For mother kernels 
LRV_mother_estimator <- function(try_b, all_autocovariances, the_kernel){
  big_T <- nrow(all_autocovariances)
  
  # Make the wieghts that correspond to the autocovariances
  # Each row corresponds to the wieghts for a specific wieght. 
  # Each column corresponds to the autocovariance lag
  # Each simulation gets each row of these scenarios. 
  W <- sapply(try_b, function(b){
    if(b == 0){
      new_weights <- c(1, rep(0, c(big_T-1))) 
    } else{
      M <- b*big_T 
      new_weights <- sapply(0:(big_T -1)/(M), the_kernel)
    }
  })
  W <- t(W)
  rownames(W) = paste("b=", round(try_b, 3), sep = "")
  colnames(W) = paste("Lag=", 0:(nrow(all_autocovariances)-1), sep = "")
  
  
  # [#, ] the try_b value 
  # [ , #]  a component of the estimated LRV 
  #    Omega_11, Omega_12,  ..., Omega_1d; ...; Omega_d1, ...Omega_dd
  omega_hats <- W %*% all_autocovariances
  rownames(omega_hats) <- paste("b=", round(try_b, 3))
  
  return(omega_hats)
}



# For Lugsail Kernels 
# Corrects the non-positive diagonal elements. 

LRV_estimator <- function(try_b, all_autocovariances, 
                          the_kernel, lugsail_parameters = list(r = 1, c= 0), 
                          mother_omega){
  d_max <- sqrt(ncol(all_autocovariances))
  big_T <- nrow(all_autocovariances)
  
  # Make the wieghts that correspond to the autocovariances
  # Each row corresponds to the weights for a specific weight. 
  # Each column corresponds to the autocovariance lag
  # Each simulation gets each row of these scenarios. 
  W <- sapply(try_b, function(b){
    if(b == 0){
      new_weights <- c(1, rep(0, c(big_T-1))) 
    } else{
      M <- b*big_T 
      new_weights <- sapply(0:(big_T -1)/(M), lugsail, the_kernel = the_kernel, 
                            lugsail_parameters = lugsail_parameters)
    }
  })
  W <- t(W)
  rownames(W) = paste("b=", round(try_b, 3), sep = "")
  colnames(W) = paste("Lag=", 0:(nrow(all_autocovariances)-1), sep = "")
  
  # [#, ] the try_b value 
  # [ , #]  a component of the estimated LRV 
  #    Omega_11, Omega_12,  ..., Omega_1d; ...; Omega_d1, ...Omega_dd
  omega_hats <- W %*% all_autocovariances
  rownames(omega_hats) <- paste("b=", round(try_b, 3))
  
  #  ---- Fix PSD Issues ----
  # diagonal indices, the columns to check for omega_hats  
  check_index <- seq(1, d_max^2, by = (d_max+1)) 
  
  # Matrix with each row corresponding to a b, and each column corresponding 
  # to a diagonal element
  # [#,  ]: try_b
  # [ , #]: diagonal element
  counts <- matrix(0, nrow = nrow(omega_hats), ncol = length(check_index))
  
  for(i in 1:length(check_index)){
    col_index <- check_index[i]
    index_fix <- omega_hats[, col_index] <= 0
    counts[index_fix ,i] <- 1  # Corrected values =1 , not-corrected = 0 
    omega_hats[index_fix, col_index] <- mother_omega[index_fix,col_index]    
  }
  
  return(omega_hats)
}




# Calculate F-Statistics --------------------------------------------------


F_stats <- function(big_T, the_means, omega_hats, d, the_Ms, null_means = rep(0, d)){
  
  # [#]: try_b values 
  F_stat_by_b <- apply(omega_hats, 1, function(one_omega_hat){
    one_omega_hat <- matrix(one_omega_hat, nrow = (d+1), ncol = (d+1))
    one_omega_hat <- one_omega_hat
    one_omega_hat <- as.matrix(nearPD(one_omega_hat)$mat)
    one_omega_hat <- (solve(the_Ms)%*%one_omega_hat%*%solve(the_Ms))[2,2] # second coef only
    inv_one_omega <- 1/one_omega_hat
    num <- (the_means- null_means)[2]
    F_stat = (num*inv_one_omega*num)*big_T
    return(F_stat)
  })
  
  
  
  return(F_stat_by_b)
}




get_kernel_F_stats <- function(try_b, the_means, d, all_autocovariances, 
                               the_kernel, q = 1, the_Ms){
  
  big_T <- nrow(all_autocovariances)
  
  # Mother Kernel 
  omega_mother <- LRV_mother_estimator(try_b, all_autocovariances, the_kernel)
  F_mother <- t(F_stats(big_T, the_means, omega_mother,  d, the_Ms))
  
  # Zero Lugsail 
  lug_para <- get_lugsail_parameters(big_T, q = q, method = "Zero")
  omega_zero <- LRV_estimator(try_b, all_autocovariances, 
                              the_kernel = the_kernel, 
                              lugsail_parameters = lug_para, 
                              omega_mother)
  F_zero <- t(F_stats(big_T, the_means, omega_zero,  d, the_Ms))
  
  
  return(list(zero = F_zero, mother = F_mother))
}




# Main --------------------------------------------------------------------



# Generates a null data set.  Calculates the autocovariance matrices. 
# Calculates the F-statistic for dimensions (1,..., d_max). 
# Need to repeat this simulation many times to then find what asymptotic CV is. 


simulate_f_stat <- function(big_T = 1000, alpha = alpha, d = 5, model_type, rho, 
                            theta = rep(0, d)){
  
  
  # Simulate the data 
  all_sim_data <- errors_and_thetas(big_T, rho, model_type, d = d, theta = theta)
  
  est_thetas <- all_sim_data[[1]] 
  orig_errors <- all_sim_data[[2]]
  the_Ms <- all_sim_data[[3]]
  
  the_means <- est_thetas      # These are the parameters estimated 
  all_sim_data <- orig_errors  # These are the errors 
  
  
  # ------- AutoCovariance Matrices  -------
  # [#, ] the lag (0, ..., big_T-1)
  # [ , #]  a component of the vectorized autocov matrix. 
  #    R11, R12,  ..., R1d; R21, R22, ..., R2d; ...; Rd1, ...Rdd
  all_autocovariances <-sapply(0:(big_T-1), R, the_sim_data = all_sim_data)
  all_autocovariances <- t(all_autocovariances)
  
  
  # ------- Autocorrelation  -------  
  # [#,  ] the lag, (0, ..., big_T-1)
  # [ , #]  each simulation, 1:num_sims 
  all_autocorrelation <- stats::acf(all_sim_data[,2], plot = F)$acf
  
  # ------- Bandwidth Rules -------  
  # [#,  ] the lag, (0, ..., big_T-1)
  # [ , #]  each simulation, 1:num_sims 
  b_mse <- mse_rule(big_T, rho = all_autocorrelation[2], q = 1)
  b_ft <- ft_rule(all_autocorrelation, big_T = big_T )
  b_opt <- opt_rule(abs(all_autocorrelation[2]), big_T = big_T, alpha = alpha, d= 1)
  try_b <- c(mse = b_mse, ft = b_ft, opt = b_opt)
  
  
  # ------- F-statistics for various b values (an 2-dimensional array) -------
  # [ #,  ] a different b value
  # [  , #] a different d
  
  # ------- BARTLETT ------- 
  F_bartlett <- get_kernel_F_stats(try_b, the_means, 4, 
                                   all_autocovariances, 
                                   bartlett, q = 1, the_Ms)
  F_mother <- F_bartlett$mother
  F_zero <- F_bartlett$zero
  
  # ------- Fixed-b Rejection Rate -------
  cv <- grid_cv_find(d = 1, alpha = alpha, the_kernel="Bartlett", 
                     lug_type = "Zero", new_b= try_b)
  rr_fixed_zero <- F_zero > cv
  
  cv <- grid_cv_find(d = 1, alpha = alpha, the_kernel="Bartlett", 
                     lug_type = "Mother", new_b= try_b)
  rr_fixed_mother <- F_mother > cv
  
  # ------- Small-b Rejection Rate -------
  cv <- qchisq((1-alpha), 1)
  rr_small_zero <- F_zero > cv
  rr_small_mother <- F_mother > cv
  
  rr <- c( fixed_mother = rr_fixed_mother,fixed_zero =rr_fixed_zero , 
           small_mother = rr_small_mother, small_zero= rr_small_zero )
  
  names(rr) <- paste(rep(c("mother", "zero"), each = 3, 2),
                     rep(c("fixed", "small"), each = 6), 
                     rep(c("mse", "ft", "opt"),4), sep = "_")
  
  return(rr)
}



# Simulation Parameters (Adjust me) ---------------------------------------


# How many replicates:  10,0000
num_replicates <- 1000

# Sample size of each replicate
big_T <- 200

# Dimensions
d <- 4

rho <- 0
model_type = "AR1_HOMO"
alpha = .05
the_kernel = bartlett 
q = 1

# Run Simulation  ---------------------------------------------------------
set.seed(26)

repeated_sims <- function(num_replicates, big_T, alpha, d, rho, model_type){
  test_stats <- replicate(num_replicates, 
                          simulate_f_stat(big_T = big_T, alpha = alpha, 
                                          d = d, rho = rho,
                                          model_type = model_type))
  
  
  
  results <- data.frame(do.call(rbind, strsplit(rownames(test_stats), "_")))
  colnames(results) <- c("lug_type", "cv", "rule")
  results$rate <- colMeans(t(test_stats))
  results$num_sims <- num_replicates
  results$big_T <- big_T
  results$rho <- rho 
  results$model <- model_type
  
  
  print(paste(model_type, big_T, rho))
  #zero_results <- results[results$lug_type == "zero", ]
  print(results[order(results$lug_type,results$rate, decreasing = T), ])
  return(results)
}


all_results <- data.frame(NA, NA, NA, NA, NA, NA, NA, NA)
colnames(all_results) <- c("lug_type", "cv", "rule", "rate", "big_T", "rho", "model", "num_sims")

for(big_T in c(200, 500, 1000)){
  for(rho in c(0.9, 0.75, 0.5, .25, 0)){
    results_M1 <- repeated_sims(num_replicates, big_T, alpha, d, rho, "AR1_HOMO")
    results_M2 <- repeated_sims(num_replicates, big_T, alpha, d, rho, "AR1_HET")
    results_M3 <- repeated_sims(num_replicates, big_T, alpha, d, rho, "ARMA_G")
    results_M4 <- repeated_sims(num_replicates, big_T, alpha, d, rho, "ARMA_L")
    all_results <- rbind(all_results, results_M1, results_M2, results_M3, results_M4) 
  }
  write.csv(all_results, file = "compare_procedures_zero_rule2.csv", row.names = F)
}

all_results <- all_results[-1, ]
write.csv(all_results, file = "compare_procedures_zero_rule2.csv", row.names = F)
