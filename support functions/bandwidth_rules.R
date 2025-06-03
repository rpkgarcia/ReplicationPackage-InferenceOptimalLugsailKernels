
# Andrews (1991) Rule as described by Lazarus2018
mse_rule <- function(big_T, rho = .25, q = 1, d = 1){
  #the_b <- 0.75*big_T^(-2*q/(2*q+1))
  w_1 <- (2*rho/(1-rho)^2)
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
  cv <- qchisq((1-alpha), d)
  num2 <- 0 
  w_q <- 2*rho/(1-rho^2)
  
  if (lug_type == "Over"){
    b_init <- opt_rule(rho, big_T, alpha, d, "Zero")
    g_q <- -1
    num2 <-  dchisq(cv, d)*cv*g_q*w_q/(b_init*big_T)
  }
  
  num1 <-  alpha/(big_T*log(rho))
  #num1 <-  .5/(big_T*log(rho))
  den1 <- dchisq(cv, d)*cv*2*rho^2/(1+ rho)
  the_b <- log(-(num1 + num2)/den1)/(big_T*log(rho))
  the_b <- max(0, the_b)
  if(is.na(the_b)){the_b <- 0 }
  
  return(the_b)
}

