# ------------------------------------------------------------------------
  # Kernel functions --------------------------------------------------------
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
