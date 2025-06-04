
library(Matrix)

# Load Functions ----------------------------------------------------------


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




# Load Functions ----------------------------------------------------------



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


# LRV Estimator -----------------------------------------------------------


# For mother kernels 
LRV_mother_estimator <- function(try_b, all_autocovariances, the_kernel){
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
  
  return(list(omega_hats, counts))
}



# Calculate F-Statistics --------------------------------------------------


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

F_stats <- function(the_means, omega_hats, d_vec, 
                    null_means = rep(0, max(d_vec))){
  
  # [#,  ]: d_vec value 
  # [ , #]: try_b values 
  F_stat_by_d_b <- apply(omega_hats, 1, function(one_omega_hat){
    F_stat_by_d <- sapply(d_vec, function(the_d){
      one_omega_hat <- matrix(one_omega_hat, nrow = sqrt(length(one_omega_hat)))
      one_omega_hat <- one_omega_hat[1:the_d, 1:the_d]
      one_omega_hat <- as.matrix(nearPD(one_omega_hat)$mat)
      inv_one_omega <- solve(one_omega_hat) 
      num <- (the_means- null_means)[1:the_d]
      F_stat = (num%*% inv_one_omega  %*%num)/the_d
    })
    names(F_stat_by_d) <- d_vec 
    return(F_stat_by_d)
  })
 
  
  return(F_stat_by_d_b)
}



# Generate all F_stats by Kernel Family -----------------------------------


get_kernel_F_stats <- function(try_b, the_means, d_vec, all_autocovariances, 
                               the_kernel, q = 1){
  
  # Mother Kernel 
  omega_mother <- LRV_mother_estimator(try_b, all_autocovariances, the_kernel)
  F_mother <- t(F_stats(the_means, omega_mother, d_vec))
  
  # Zero Lugsail 
  lug_para <- get_lugsail_parameters(big_T, q = q, method = "Zero")
  omega_zero <- LRV_estimator(try_b, all_autocovariances, 
                              the_kernel = the_kernel, 
                              lugsail_parameters = lug_para, 
                              omega_mother)
  F_zero <- t(F_stats(the_means, omega_zero[[1]], d_vec))
  
  # # Adaptive Lugsail 
  # big_T <- nrow(all_autocovariances)
  # init_b <-  0.75*big_T^(-2*q/(2*q+1))
  # lug_para <- get_lugsail_parameters(big_T, q = q, method = "Adaptive", 
  #                                    b = init_b)
  # omega_adapt <- LRV_estimator(try_b, all_autocovariances, 
  #                              the_kernel= the_kernel, 
  #                              lugsail_parameters = lug_para, 
  #                              omega_mother)
  # F_adapt <- t(F_stats(the_means, omega_adapt[[1]], d_vec))
  # 
  # # Over Lugsail 
  # lug_para <- get_lugsail_parameters(big_T, q = q, method = "Over")
  # omega_over <- LRV_estimator(try_b, all_autocovariances, 
  #                             the_kernel= the_kernel,
  #                             lugsail_parameters = lug_para, 
  #                             omega_mother)
  # F_over <- t(F_stats(the_means, omega_over[[1]], d_vec))
  
  return(list(F_stats = list(mother = F_mother, zero = F_zero), 
                             #adapt = F_adapt, over = F_over), 
              psd_counts = data.frame(zero = omega_zero[[2]])))#, 
                                      #adapt = omega_adapt[[2]], 
                                      #over = omega_over[[2]])))
}


# Main --------------------------------------------------------------------



# Generates a null data set.  Calculates the autocovariance matrices. 
# Calculates the F-statistic for dimensions (1,..., d_max). 
# Need to repeat this simulation many times to then find what asymptotic CV is. 


simulate_f_stat <- function(big_T = 1000, d_vec = c(1, 5, 10)){
  d_max <- max(d_vec)
  
  # Simulate the data 
  sim_data <- matrix(rnorm(d_max*big_T), nrow = big_T, ncol = d_max)
  the_means <- colMeans(sim_data)
  sim_data <- apply(sim_data, 1, function(row) row - the_means)
  sim_data <- t(sim_data)
  
  # ------- AutoCovariance Matrices  -------
  # [#, ] the lag (0, ..., big_T-1)
  # [ , #]  a component of the vectorized autocov matrix. 
  #    R11, R12,  ..., R1d; R21, R22, ..., R2d; ...; Rd1, ...Rdd
  all_autocovariances <-sapply(0:(big_T-1), R, the_sim_data = sim_data)
  all_autocovariances <- t(all_autocovariances)

  
  
  # ------- F-statistics for various b values (an 2-dimensional array) -------
  # [ #,  ] a different b value
  # [  , #] a different d
  
  # ------- BARTLETT ------- 
  F_bartlett <- get_kernel_F_stats(try_b, the_means, d_vec, all_autocovariances, 
                                   bartlett, q = 1)
  psd_bartlett <- F_bartlett[[2]]
  F_bartlett <- F_bartlett[[1]]
  
  # # ------- PARZEN -------
  # F_parzen <- get_kernel_F_stats(try_b, the_means, d_vec, all_autocovariances,
  #                                parzen, q = 2)
  # psd_parzen <- F_parzen[[2]]
  # F_parzen <- F_parzen[[1]]
  # 
  # # -------  QS -------
  # F_qs <- get_kernel_F_stats(try_b, the_means, d_vec, all_autocovariances,
  #                            qs, q = 2)
  # psd_qs <- F_qs[[2]]
  # F_qs <- F_qs[[1]]
  
  
  
  # ------- Store all F Stats -------
  # [ #,  ,  ] a different b value (try_b)
  # [  , #,  ] a different d value (d_vec)
  # [  ,  , #] a the_kernel 
  F_stats_all <- array(c(F_bartlett[["mother"]], F_bartlett[["zero"]]), 
                         #F_bartlett[["adapt"]], F_bartlett[["over"]]), 
                         # F_parzen[["mother"]], F_parzen[["zero"]], 
                         # F_parzen[["adapt"]], F_parzen[["over"]], 
                         # F_qs[["mother"]], F_qs[["zero"]], 
                         # F_qs[["adapt"]], F_qs[["over"]]),
                       dim=c(length(try_b), length(d_vec), 2))
  
  
  dimnames(F_stats_all)[[3]] <-  c("Bartlett_Mother", "Bartlett_Zero")#,
                                   # "Bartlett_Over") #  "Bartlett_Adapt",
                                   # "Parzen_Mother", "Parzen_Zero",
                                   #  "Parzen_Adapt", "Parzen_Over", 
                                   #  "QS_Mother", "QS_Zero",
                                   #  "QS_Adapt", "QS_Over")
  
  # ------- Store PSD Statistics -------
  psd_counts <- list(bartlett = psd_bartlett) 
                     # parzen = psd_parzen, 
                     # qs = psd_qs)
  
  return(list(F_stats_all = F_stats_all, psd_counts = psd_counts))
}


# Simulation Parameters (Adjust me) ---------------------------------------


#  What bs to use (should be 0.005 to 0.99)
try_b <-  seq(0.005, .5, by = 0.005)

# How many replicates
# KV005 used 50,0000
num_replicates <- 50000

# Sample size of each replicate
big_T = 1000

# Maximum number of dimensions
# Should do 12
d_vec <- c(2:4)


# -------------------------------------------------------------------------
# Run Simulation  ---------------------------------------------------------
# -------------------------------------------------------------------------
set.seed(26)
test_stats <- replicate(num_replicates, 
                        simulate_f_stat(big_T = big_T, d_vec = d_vec))

psd_counts <- test_stats[2,]
#all_F_stats <- test_stats[1, ]

#  ---- Save count info ------ 
get_counts <- function(the_mother_kernel){
  counts_matrix<- matrix(0, nrow = length(try_b), ncol = max(d_vec)*3) 
  lug_type <- colnames(psd_counts[[1]][[the_mother_kernel]])
  
  for(index in 1:num_replicates){
      counts_matrix <- counts_matrix + psd_counts[[index]][[the_mother_kernel]]
  } 
  counts_zero <- counts_matrix[,grepl("zero", lug_type)]
  #counts_adapt <- counts_matrix[,grepl("adapt", lug_type)]
  #counts_over <- counts_matrix[,grepl("over", lug_type)]
  psd_counts <- list(zero = counts_zero) 
   #                  adapt = counts_adapt, 
                    # over =counts_over) 
  return(psd_counts)
}

all_psd_counts <- list(bartlett = get_counts("bartlett")) 
                       # parzen = get_counts("parzen"), 
                       # qs = get_counts("qs"))


# ----- Save test statistic info ------ 
all_F_stats <- c()
for(i in seq(1, num_replicates*2, by = 2)){
  all_F_stats  <- c(all_F_stats, test_stats[[i]])
}

all_F_stats <- array(c(all_F_stats),
                     dim = c(length(try_b), length(d_vec), 2, num_replicates))


# [ #,  ,  ,  ] : rows, the different b values [1:length(try_b)]
# [  , #,  ,  ] : columns, the different d values [1:max_d]
# [  ,  , #,  ] : the kernel type 
# [  ,  ,  , #] : each simulation [1:num_replicates]

dimnames(all_F_stats)[[1]] <- try_b
dimnames(all_F_stats)[[2]] <- d_vec
dimnames(all_F_stats)[[3]] <-  c("Bartlett_Mother", "Bartlett_Zero")
                                 #"Bartlett_Adapt", "Bartlett_Over") 
                                 # "Parzen_Mother", "Parzen_Zero",
                                 # "Parzen_Adapt", "Parzen_Over", 
                                 # "QS_Mother", "QS_Zero",
                                 # "QS_Adapt", "QS_Over")
dimnames(all_F_stats)[[4]] <- paste("sim", 1:num_replicates, sep="")

# -------------------------------------------------------------------------
# Create CV Tables --------------------------------------------------------
# -------------------------------------------------------------------------

save_cv <- function(the_kernel){
  F_stats <- all_F_stats[ , , the_kernel, ]
  
  t90 <- apply(F_stats, 1:2, quantile, probs = 0.90)
  write.csv(t90, paste(the_kernel, "_10.csv", sep = ""))
  
  t95 <- apply(F_stats, 1:2, quantile, probs = 0.95)
  write.csv(t95, paste(the_kernel, "_05.csv", sep = ""))
  
  t97.5 <- apply(F_stats, 1:2, quantile, probs = 0.975)
  write.csv(t97.5, paste(the_kernel, "_025.csv", sep = ""))
  
  t99 <- apply(F_stats, 1:2, quantile, probs = 0.99)
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

save_distr <- function(the_kernel, d){
  F_stats <- t(all_F_stats[ , d, the_kernel, ])
  write.csv(F_stats, paste(the_kernel, "_d", d_vec[d], ".csv", sep = ""))
}

#  this is how I get all the critical values 
setwd("Fixed_b_distribution")
sapply(1:dim(all_F_stats)[2], function(d){
  sapply(dimnames(all_F_stats)[[3]], save_distr, d = d)
})
setwd("..")



# -------------------------------------------------------------------------
# PSD Plots ---------------------------------------------
# -------------------------------------------------------------------------

# Saves a data frame for a specific the_kernel and dimension. 
# Keeps the F-statistics so I can use them to make distribution plots later. 

save_psd_info <- function(the_kernel){
  sapply(names(all_psd_counts[[the_kernel]]), function(lug_type){
    psd_info <- all_psd_counts[[the_kernel]][[lug_type]]/num_replicates
    write.csv(psd_info, paste(the_kernel,"_",lug_type, ".csv", sep = ""))
  })
}

#  this is how I get psd rate_plots 
setwd("psd_issues")
sapply(names(all_psd_counts), save_psd_info)
setwd("..")

