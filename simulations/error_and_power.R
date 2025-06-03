library(Matrix)
library(MASS)



# Get CV ------------------------------------------------------------------

grid_cv_find <- function(d, alpha, the_kernel, lug_type, new_b){
  # Make alpha a character
  alpha <- alpha*100
  if(alpha == 10){
    alpha <- "10"
  } else {
    alpha <- paste("0", alpha, sep = "")
  }
  
  file <- paste(the_kernel, lug_type, alpha, sep = "_")
  grid <- read.csv(paste("../../fixed_b_asymptotics/",file, "_Master.csv", sep = ""))
  
  cv_by_b <- sapply(new_b, function(b){
    return(grid[which.min(abs(b-grid[,1])),d+1])
  })
  
  return(cv_by_b)
}


# Simulate a AR(1) Model -------------------------------------------------


library(MASS)

# Call outside of function for computation time. 
exponent <- abs(matrix(1:10 - 1, nrow = 10, ncol = 10, byrow = TRUE)-(1:10 - 1))
Sigma <- 0.9^exponent

AR_X <- function(big_T, rho, d){
  X <- matrix(0, nrow = big_T, ncol = d)
  
  # Sigma[1:d, 1:d]
  for(i in 2:big_T){
    X[i, ] <- rho*X[(i-1), ] + mvrnorm(n = 1, rep(0, d), diag(d)) 
  }
  return(X)
}

VAR1_SIMPLE <- function(big_T, rho, d, mean_vec = 0 ){
  u <- AR_X(big_T, rho, d)
  Y <- mean_vec + u
  return(list(Y = Y))
}


# Practice calling 
errors_and_thetas <- function(big_T,  d = 1, rho = 0, mean_vec = 0){
  the_data <- VAR1_SIMPLE(big_T, rho, d, mean_vec)
  X <- rep(1, big_T)
  
  coefs <- (solve(t(X)%*%X)%*% t(X)%*%the_data$Y)  # All coefs 
  M <- (t(X)%*%X/big_T) 
  errors <- X*(the_data$Y - X%*%coefs)
  errors <- errors - apply(errors, 2, mean)  # Center them immediately 
  
  return(list(coefs= coefs, errors = errors, M = M))
}



# Load Functions ----------------------------------------------------------
# source("../../support_functions/kernels.R")

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

# source("../../support_functions/est_autocov.R")

R_multi <- function(h, the_sim_data){
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



R_uni <- function(){
  
}

R <- function(h, the_sim_data){
  return(R_multi(h, the_sim_data))
}

# source("../../support_functions/error_rates.R")
# source("../../support_functions/bandwidth_rules.R")
# Andrews (1991) Rule as described by Lazarus2018
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
  num1 <-  alpha^(1/(2*d))/(big_T*log(rho))
  den1 <- dchisq(cv, d)*cv*2*rho^2/(1+ rho)
  the_b <- log(-(num1 + num2)/den1)/(big_T*log(rho))
  the_b <- max(0, the_b)
  if(is.na(the_b)){the_b <- 0 }
  
  return(the_b)
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
    if(b == 0 ){
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
  # Each row corresponds to the wieghts for a specific wieght. 
  # Each column corresponds to the autocovariance lag
  # Each simulation gets each row of these scenarios. 
  W <- sapply(try_b, function(b){
    if(b == 0 ){
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


F_stats <- function(the_means, omega_hats, d, the_Ms, null_means = rep(0, d), big_T){
  
  # [ #]: try_b values 
  F_stat_by_b <- apply(omega_hats, 1, function(one_omega_hat){
    one_omega_hat <- matrix(one_omega_hat, nrow = (d), ncol = (d))
    one_omega_hat <- one_omega_hat
    one_omega_hat <- as.matrix(nearPD(one_omega_hat)$mat)
    inv_one_omega <- solve(one_omega_hat)
    num <- (the_means- null_means)
    F_stat = (num%*%inv_one_omega%*%t(num))*big_T/d
    return(F_stat)
  })
  
  
  
  return(F_stat_by_b)
}




get_kernel_F_stats <- function(try_b, the_means, d, all_autocovariances, 
                               the_kernel, q = 1, the_Ms){
  
  big_T <- nrow(all_autocovariances)
  
  # Mother Kernel 
  omega_mother <- LRV_mother_estimator(try_b, all_autocovariances, the_kernel)
  F_mother <- t(F_stats(the_means, omega_mother,  d, the_Ms, null_means =rep(0, d), 
                        big_T = big_T))
  
  # Zero Lugsail 
  lug_para <- get_lugsail_parameters(big_T, q = q, method = "Zero")
  omega_zero <- LRV_estimator(try_b, all_autocovariances, 
                              the_kernel = the_kernel, 
                              lugsail_parameters = lug_para, 
                              omega_mother)
  F_zero <- t(F_stats(the_means, omega_zero,  d, the_Ms, null_means =rep(0, d), 
                      big_T = big_T))
  
  
  return(list(zero = F_zero, mother = F_mother))
}




# Main --------------------------------------------------------------------



# Generates a null data set.  Calculates the autocovariance matrices. 
# Calculates the F-statistic for dimensions (1,..., d_max). 
# Need to repeat this simulation many times to then find what asymptotic CV is. 


simulate_f_stat <- function(big_T=1000, alpha=0.05, d=5, rho=0, mean_vec=rep(0, d), try_b){
  
  # Simulate the data 
  all_sim_data <- errors_and_thetas(big_T, d=d, rho=rho, mean_vec=mean_vec)
  
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
  
  
  # ------- F-statistics for various b values (an 2-dimensional array) -------
  # [ #,  ] a different b value
  # [  , #] a different d
  
  # ------- BARTLETT ------- 
  F_bartlett <- get_kernel_F_stats(try_b, the_means, d, 
                                   all_autocovariances, 
                                   bartlett, q = 1, the_Ms)
  F_mother <- F_bartlett$mother
  F_zero <- F_bartlett$zero
  
  
  return(c(mother = F_mother, zero = F_zero))
}




# Simulation Parameters (Adjust me) ---------------------------------------


# How many replicates:  1000
num_replicates <-1000

# Sample size of each replicate
big_T = 500

model_type = "AR1_HOMO"
alpha = .05



#  What bs to use 
try_b <-  c(seq(0, .20, by= 0.01),
            #seq(0.11, .20, by = 0.01),
            seq(0.25, .5, by =0.05), 
            mse_rule(big_T = big_T, rho = 0.75), 
            opt_rule(rho = 0.75, big_T = big_T, alpha = 0.05, d= 2),
            opt_rule(rho = 0.75, big_T =big_T, alpha = 0.05, d= 3),
            opt_rule(rho = 0.75, big_T = big_T, alpha = 0.05, d= 4))
try_b <- sort(try_b)

# Generate data -----------------------------------------------------------


# d<-2
# rho <- 0 
# recommended_delta_sq <- c(2.4306,3.2463,3.8148, 4.2728, 4.6658,5.0149, 5.3321,
#                           5.6245, 5.8973, 6.1538)
# #delta_sq <- recommended_delta_sq[d]
# TRUE_LRV <- 1/(1 - rho^2)
# under_alt <- sqrt(recommended_delta_sq[d]/(big_T*TRUE_LRV))  <- this was wrong!!!! 



# Save Information --------------------------------------------------------


run_and_save <- function(num_replicates, big_T, alpha, d, try_b, rho, mean_vec="Null"){
  recommended_delta_sq <- c(2.4306, 3.2463,3.8148, 4.2728, 4.6658,5.0149, 5.3321, 5.6245, 5.8973, 6.1538)
  
  if(mean_vec == "Null"){
    mean_vec <- 0 
  } else if (mean_vec == "Alt") {
    mean_vec <- sqrt(recommended_delta_sq[d]/(1 - rho^2)/(big_T))
  }
  
  test_stats <- replicate(num_replicates, 
                          simulate_f_stat(big_T = big_T,alpha = alpha, 
                                          d = d, rho = rho, mean_vec = mean_vec, 
                                          try_b = try_b))
  
  F_mother <- test_stats[grepl("mother", rownames(test_stats)),]
  F_zero <- test_stats[grepl("zero", rownames(test_stats)),]
  
  num_replicates <- ncol(test_stats)
  rownames(F_mother) <- try_b
  rownames(F_zero) <- try_b
  
  cv_zero <-  grid_cv_find(d, alpha, "Bartlett", "Zero", try_b)
  rr_zero <-rowSums(apply(F_zero, 2, function(sim) sim>cv_zero))/num_replicates
  
  cv_mother <- grid_cv_find(d, alpha, "Bartlett", "Mother", try_b)
  rr_mother <-rowSums(apply(F_mother, 2, function(sim) sim>cv_mother))/num_replicates
  
  cv_small <- qchisq(.95, d)/d
  rr_mother_small <- rowSums(apply(F_mother, 2, function(sim) sim>cv_small))/num_replicates
  rr_zero_small <- rowSums(apply(F_zero, 2, function(sim) sim>cv_small))/num_replicates
  
  # Save Everything 
  lug_types <- rep(c("Zero", "Mother", "Zero", "Mother"), each = length(try_b)) 
  cv_type <- rep(c("fixed", "small"), each = length(try_b)*2)
  rr <- c(rr_zero, rr_mother, rr_zero_small, rr_mother_small)
  results <- data.frame(lug_type = lug_types, 
                        b = try_b, 
                        ncp = mean_vec[1], 
                        cv_type = cv_type, 
                        rr = rr, 
                        rho = rho,
                        d = d,
                        big_T= big_T, 
                        num_sims = num_replicates, 
                        use = "error_rates")
  print(paste("rho", rho, "; d", d, "; mean_vec", mean_vec))
  return(results)
}

# Under the NULL Hypothesis
test_stats_d2_r0 <- run_and_save(num_replicates = num_replicates, big_T = big_T, alpha = 0.05, d = 2, try_b = try_b, rho = 0)
test_stats_d2_r75 <- run_and_save(num_replicates = num_replicates, big_T = big_T, alpha = 0.05, d = 2, try_b = try_b, rho = 0.75)

test_stats_d3_r0 <- run_and_save(num_replicates = num_replicates, big_T = big_T, alpha = 0.05, d = 3, try_b = try_b, rho = 0)
test_stats_d3_r75 <- run_and_save(num_replicates = num_replicates, big_T = big_T, alpha = 0.05, d = 3, try_b = try_b, rho = .75)

test_stats_d4_r0 <- run_and_save(num_replicates = num_replicates, big_T = big_T, alpha = 0.05, d = 4, try_b = try_b, rho = 0)
test_stats_d4_r75 <- run_and_save(num_replicates = num_replicates, big_T = big_T, alpha = 0.05, d = 4, try_b = try_b, rho = .75)



# Under the ALT Hypothesis
recommended_delta_sq <- c(2.4306, 3.2463,3.8148, 4.2728, 4.6658,5.0149, 5.3321, 5.6245, 5.8973, 6.1538)
#mean_vec = sqrt(recommended_delta_sq[d]/(1 - rho^2)/(big_T))

test_stats_alt_d2_r0 <- run_and_save(num_replicates = num_replicates, big_T = big_T, alpha = 0.05, d = 2, try_b = try_b, rho = 0, mean_vec = "Alt")
test_stats_alt_d2_r75 <- run_and_save(num_replicates = num_replicates, big_T = big_T, alpha = 0.05, d = 2, try_b = try_b, rho = 0.75, mean_vec = "Alt")


test_stats_alt_d3_r0 <- run_and_save(num_replicates = num_replicates, big_T = big_T, alpha = 0.05, d = 3, try_b = try_b, rho = 0, mean_vec = "Alt")
test_stats_alt_d3_r75 <- run_and_save(num_replicates = num_replicates, big_T = big_T, alpha = 0.05, d = 3, try_b = try_b, rho = .75, mean_vec = "Alt")


test_stats_alt_d4_r0 <- run_and_save(num_replicates =num_replicates, big_T = big_T, alpha = 0.05, d = 4, try_b = try_b, rho = 0, mean_vec = "Alt")
test_stats_alt_d4_r75 <- run_and_save(num_replicates = num_replicates, big_T = big_T, alpha = 0.05, d = 4, try_b = try_b, rho = .75, mean_vec = "Alt")



# Power -------------------------------------------------------------------
power_rates <- function(test_stats, d, rho, try_b){
  F_mother <- c(test_stats[grepl("mother", rownames(test_stats)),])
  F_zero <- c(test_stats[grepl("zero", rownames(test_stats)),])
  
  cv_zero <-  grid_cv_find(d, alpha, "Bartlett", "Zero", try_b["Zero"])
  rr_zero_fixed <-  sum(F_zero>cv_zero)/length(F_zero)
  
  cv_mother <- grid_cv_find(d, alpha, "Bartlett", "Mother", try_b["Mother"])
  rr_mother_fixed <- sum(F_mother>cv_mother)/length(F_mother)
  
  cv_small <- qchisq(.95, d)/d
  rr_mother_small <- sum(F_mother>cv_small)/length(F_mother)
  rr_zero_small <- sum(F_zero>cv_small)/length(F_zero)
  
  rr <- c("fixed_zero" = rr_zero_fixed, 
          "fixed_mother" = rr_mother_fixed,
          "small_zero" = rr_zero_small, 
          "small_mother" = rr_mother_small)
  return(rr)
}



# Power by Sample Size ----------------------------------------------------

run_and_save_power <- function(num_replicates= 50, big_T, rho= 0, d = 2){
  try_delta_sq <- seq(0, 12, by = 3)
  try_b = c("Zero"= opt_rule(rho, big_T, 0.05 , d, "Zero"), 
            "Mother" = mse_rule(big_T = big_T, rho = rho))
  
  #delta_sq <- recommended_delta_sq[d]
  TRUE_LRV <- 1/(1 - rho^2)
  
  power_test_stats <-data.frame(NA, NA, NA, NA, NA, NA, NA, NA, NA,  NA)
  colnames(power_test_stats) <- c("lug_type", "b", "ncp", "cv_type", "rr","rho","d","big_T", "num_sims",  "use")
  
  for(delta_sq in try_delta_sq){
    print(paste("delta_sq=",delta_sq, "; b_opt=", try_b[1], "; b_mse=",try_b[2], "; rho=", rho, "; d=", d, sep =""))
    under_alt <- sqrt(delta_sq*TRUE_LRV/big_T)
    power_F_stats <- replicate(num_replicates, 
                               simulate_f_stat(big_T = big_T,alpha = alpha, 
                                               d = d, rho = rho, mean_vec = under_alt, try_b = try_b))
    power_F_stats <- power_F_stats[c(2, 3),]
    obs_rr <- power_rates(power_F_stats, d, rho, try_b)
    
    lug_types <- c("Zero", "Mother", "Zero", "Mother")
    cv_types <- c("fixed", "fixed", "small", "small")
    
    values <- data.frame(lug_type = lug_types, 
                         b =  c(try_b, try_b), 
                         ncp = delta_sq, 
                         cv_type = cv_types, 
                         rr = obs_rr, 
                         rho = rho, 
                         d = d, 
                         big_T = big_T, 
                         num_sims = num_replicates, 
                         use = "power")
    
    power_test_stats <- rbind(power_test_stats, values)
  }
  
  power_test_stats <- power_test_stats[-1,]
  return(power_test_stats)
}

set.seed(26)
power_test_stats_d2_rho0 <- run_and_save_power(num_replicates= num_replicates, big_T = big_T, rho = 0, d = 2)
power_test_stats_d2_rho75 <- run_and_save_power(num_replicates=num_replicates,  big_T = big_T, rho = 0.75, d= 2)

power_test_stats_d3_rho0 <- run_and_save_power(num_replicates= num_replicates, big_T = big_T, rho = 0, d = 3)
power_test_stats_d3_rho75 <- run_and_save_power(num_replicates= num_replicates, big_T = big_T, rho = 0.75, d= 3)

power_test_stats_d4_rho0 <- run_and_save_power(num_replicates=num_replicates, big_T = big_T, rho = 0, d = 4)
power_test_stats_d4_rho75 <- run_and_save_power(num_replicates=num_replicates, big_T = big_T, rho = 0.75, d= 4)


# Save everything 
keep_all <- rbind(test_stats_d2_r0, test_stats_d2_r0,
                  test_stats_d3_r0, test_stats_d3_r0,
                  test_stats_d4_r0, test_stats_d4_r0,
                  test_stats_alt_d2_r0, test_stats_alt_d2_r0,
                  test_stats_alt_d3_r0, test_stats_alt_d3_r0,
                  test_stats_alt_d4_r0, test_stats_alt_d4_r0,
                  power_test_stats_d2_rho0 , power_test_stats_d2_rho75, 
                  power_test_stats_d3_rho0 , power_test_stats_d3_rho75, 
                  power_test_stats_d4_rho0 , power_test_stats_d4_rho75)


write.csv(keep_all, paste("errors_power_T", big_T, ".csv", sep = ""), row.names = F) 

# 1 x 3 plots -------------------------------------------------------------



plot_type1<- function(test_stats, d, rho, big_T = 500, the_xlim = c(0, .5), the_ylim = c(0, .5)){
  
  if(rho == 0 & d == 2){
    the_ylab <- expression(paste("d = 2; ", rho,"= 0"))
  } else if(d== 2 & rho == 0.75){
    the_ylab <- expression(paste("d = 2; ", rho,"= 0.75"))
  } else if(d== 3 & rho == 0){
    the_ylab <- expression(paste("d = 3; ", rho,"= 0"))
  } else if(d == 3 & rho == 0.75){
    the_ylab <- expression(paste("d = 3; ", rho,"= 0.75"))
  }else if(d== 4 & rho == 0){
    the_ylab <- expression(paste("d = 4; ", rho,"= 0"))
  } else if(d == 4 & rho == 0.75){
    the_ylab <- expression(paste("d = 4; ", rho, "= 0.75"))
  }
  
  rr_mother_small <- test_stats$rr[test_stats$lug_type == "Mother" & test_stats$cv_type == "small"]
  rr_mother_fixed <- test_stats$rr[test_stats$lug_type == "Mother" & test_stats$cv_type == "fixed"]
  rr_zero_small <- test_stats$rr[test_stats$lug_type == "Zero" & test_stats$cv_type == "small"]
  rr_zero_fixed <- test_stats$rr[test_stats$lug_type == "Zero" & test_stats$cv_type == "fixed"]
  
  
  plot(try_b, rr_mother_fixed, xlim = the_xlim, ylim = the_ylim, type = "l", lty = 1, lwd = 2,
       col = 2, ylab = the_ylab, main = paste("Type 1 Error"), xlab = "b")
  abline(h = 0.05, col = "grey")
  lines(try_b, rr_zero_fixed, col ="blue", lty = 1, lwd = 2)
  lines(try_b, rr_zero_small, col = "blue", lty = 2, lwd = 2)
  lines(try_b, rr_mother_small, col = 2, lty = 2, lwd = 2)
  
  if(rho >= 0){
    b_zero_index <- which.min(abs(try_b - opt_rule(rho, big_T, .05, d)))
    points(opt_rule(rho, big_T, .05, d), rr_zero_fixed[b_zero_index], pch = 21, bg = "white",col = "blue", cex = 2, lwd = 2)
    
    b_mse_index <- which.min(abs(try_b -mse_rule(big_T, rho)))
    points(mse_rule(big_T, rho), rr_mother_small[b_mse_index], pch = 24, col = 2, bg = "white", cex = 2, lwd = 2) 
  }
}

plot_type2<- function(test_stats, d, rho, big_T = 500){
  
  rr_mother_small <- 1-test_stats$rr[test_stats$lug_type == "Mother" & test_stats$cv_type == "small"]
  rr_mother_fixed <- 1-test_stats$rr[test_stats$lug_type == "Mother" & test_stats$cv_type == "fixed"]
  rr_zero_small <- 1-test_stats$rr[test_stats$lug_type == "Zero" & test_stats$cv_type == "small"]
  rr_zero_fixed <- 1-test_stats$rr[test_stats$lug_type == "Zero" & test_stats$cv_type == "fixed"]
  
  # Same for each situation
  try_b <- test_stats$b[test_stats$lug_type == "Zero" & test_stats$cv_type == "fixed"]
  
  plot(try_b, rr_mother_fixed, xlim = c(0, 0.5),ylim = c(0, 1), type = "l", lty = 1, lwd = 2,
       col = 2, main = "Type 2 Error", ylab = "", xlab = "b")
  lines(try_b, rr_zero_fixed, col = "blue", lty = 1, lwd = 2)
  lines(try_b, rr_zero_small, col = "blue", lty = 2, lwd = 2)
  lines(try_b, rr_mother_small, col = 2, lty = 2, lwd = 2)
  
  if(rho >=0){
    b_zero_index <- which.min(abs(try_b - opt_rule(rho, big_T, .05, d)))
    points(opt_rule(rho, big_T, .05, d), rr_zero_fixed[b_zero_index],pch = 21, bg = "white",col = "blue", cex = 2, lwd = 2)
    
    b_mse_index <- which.min(abs(try_b -mse_rule(big_T, rho)))
    points(mse_rule(big_T, rho), rr_mother_small[b_mse_index],pch = 24, bg = "white",col = 2, cex = 2, lwd=2) 
  }
}

plot_power<- function(test_stats, rho, big_T, alpha, d, do_legend = F){
  b_opt <- opt_rule(rho = rho, big_T = big_T, alpha = alpha, d=d)
  b_mse <- mse_rule(big_T= big_T, rho = rho)
  
  rr_mother_small <- test_stats$rr[test_stats$lug_type == "Mother" & test_stats$cv_type == "small" & test_stats$b == b_mse]
  rr_mother_fixed <- test_stats$rr[test_stats$lug_type == "Mother" & test_stats$cv_type == "fixed" & test_stats$b == b_mse]
  rr_zero_small <- test_stats$rr[test_stats$lug_type == "Zero" & test_stats$cv_type == "small"& test_stats$b == b_opt]
  rr_zero_fixed <- test_stats$rr[test_stats$lug_type == "Zero" & test_stats$cv_type == "fixed"& test_stats$b == b_opt]
  
  # Same for all 
  ncp <- test_stats$ncp[test_stats$lug_type == "Zero" & test_stats$cv_type == "fixed"& test_stats$b == b_opt]
  
  plot(ncp, rr_mother_small, col = 2,  lty = 2,  lwd = 2, 
       type = "l", ylim = c(0, 1),  main = "Power", ylab = "", xlab = expression(delta^2))
  lines(ncp, rr_zero_fixed, col ="blue", lty = 1, lwd = 2)
  lines(ncp, rr_zero_small, col = "blue", lty = 2, lwd = 2)
  lines(ncp, rr_mother_fixed, col = 2, lty = 1, lwd = 2)
  
  if(do_legend){
    legend("topright", 
           c("Zero-Fixed", "Zero-Small", "Mother-Fixed", "Mother-Small"), 
           col = c("blue", "blue", "red", "red"), 
           lty = c(1, 2, 1, 2), 
           lwd = c(2, 2, 2, 2))
  }
}


start_plot <- function(plot_name, height = 6, width = 8){
  reso <- 300
  length <- 4
  png(paste(plot_name, ".png", sep = ""),
      units="in", 
      res = reso, 
      width = width, 
      height =height)
}





# Plot for d = 2
start_plot(paste("T", big_T, "_d", 2, sep = ""))
par(mfrow = c(2,3))
plot_type1(test_stats_d2_r0, 2, 0, big_T=big_T)
plot_type2(test_stats_alt_d2_r0, 2, 0, big_T=big_T)
plot_power(power_test_stats_d2_rho0, 0, big_T=big_T, 0.05, 2)

plot_type1(test_stats_d2_r75, 2, 0.75, big_T=big_T)
plot_type2(test_stats_alt_d2_r75, 2, 0.75, big_T=big_T)
plot_power(power_test_stats_d2_rho75, .75, big_T=big_T, 0.05, 2)
dev.off()



# Plot for d = 3
start_plot(paste("T", big_T, "_d", 3, sep = ""))
par(mfrow = c(2, 3))
plot_type1(test_stats_d3_r0, 3, 0,big_T=big_T)
plot_type2(test_stats_alt_d3_r0,3, 0,big_T=big_T)
plot_power(power_test_stats_d3_rho0, 0,big_T=big_T, 0.05, 3)

plot_type1(test_stats_d3_r75, 3, 0.75,big_T=big_T)
plot_type2(test_stats_alt_d3_r75, 3, 0.75,big_T=big_T)
plot_power(power_test_stats_d3_rho75, .75,big_T=big_T, 0.05, 3)
dev.off()



# Plot for d = 4
start_plot(paste("T", big_T, "_d", 4, sep = ""))
par(mfrow = c(2,3))
plot_type1(test_stats_d4_r0, 4, 0,big_T=big_T)
plot_type2(test_stats_alt_d4_r0, 4, 0,big_T=big_T)
plot_power(power_test_stats_d4_rho0, 0,big_T=big_T, 0.05, 4)

plot_type1(test_stats_d4_r75, 4, 0.75,big_T=big_T)
plot_type2(test_stats_alt_d4_r75, 4, 0.75,big_T=big_T)
plot_power(power_test_stats_d4_rho75, .75,big_T=big_T, 0.05, 4)
dev.off()



# For Presentation  -------------------------------------------------------

start_plot(paste("T", big_T, "_d", 2, "pp", sep = ""), height = 5)
par(mfrow = c(1,2))
plot_type1(test_stats_d2_r75, 2, 0.75, the_xlim = c(0, .25), the_ylim = c(0, .4))
plot_power(power_test_stats_d2_rho75, .75, 500, 0.05, 2, do_legend = T)
dev.off()


start_plot(paste("T", big_T, "_d", 3, "pp", sep = ""), height = 5)
par(mfrow = c(1,2))
plot_type1(test_stats_d3_r75, 3, 0.75, the_xlim = c(0, .25), the_ylim = c(0, .4))
plot_power(power_test_stats_d3_rho75, .75, 500, 0.05, 3, do_legend = T)
dev.off()

start_plot(paste("T", big_T, "_d", 4, "pp", sep = ""), height = 5)
par(mfrow = c(1,2))
plot_type1(test_stats_d4_r75, 4, 0.75, the_xlim = c(0, .25), the_ylim = c(0, .4))
plot_power(power_test_stats_d4_rho75, .75, 500, 0.05, 4,do_legend = T)
dev.off()





