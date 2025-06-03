# Error Rates

e1 <- function(big_T, rho, b, d, lug_type, alpha = 0.05){
  
  w_q <- 2*rho/(1-rho^2)  # All (but zero lugsail)
  g_q <- 0                # default zero lugsail
  c_b <- 2*rho^2*rho^(b*big_T)/(1+ rho)
  c1 = 1.5
  
  if(lug_type == "Mother"){
    g_q = 1
    c2 =0.67
    c1 = 1
  } else if(lug_type == "Over"){
    g_q = -1
    c2 <- 2.33
    c1 = 2
  }
  
  chi_sq_cv <- qchisq((1-alpha), df = d)/d
  p1 <- alpha 
  p2 <- dchisq(chi_sq_cv, df= d)*chi_sq_cv*c_b
  p3 <- dchisq(chi_sq_cv, df= d)*chi_sq_cv*g_q*w_q*(b*big_T)^(-1)
  
  
  return((p1+p2 +p3))
}


# My Result 
e2 <- function(big_T, rho, b, d, alpha = 0.05, lug_type, delta_sq, cv=qchisq((1-alpha), df = d), sun_like = T){
  
  w_q <- 2*rho/(1-rho^2)  # All (but zero lugsail)
  g_q <- 0                # default zero lugsail
  c2 <- 1.33              
  c1 <- 1.5
  
  c_b <- 2*rho^2*rho^(b*big_T)/(1+ rho)
  
  if(lug_type == "Mother"){
    g_q = 1
    c2 =0.67
    c1 <- 1
  } else if(lug_type == "Over"){
    g_q = -1
    c2 <- 2.33
    c1 <- 2
  }
  
  #cv <- grid_cv_find(d, alpha, "Bartlett", "Zero", b)
  
  p1 <- pchisq(cv, df =d, ncp = delta_sq)
  p2 <- -exp(-delta_sq/2)*dchisq(cv, df=d)*cv*.5*(cv+2-d)*c2*b
  p3 <- dchisq(cv, df = d, ncp = delta_sq)*cv*c_b
  # p4.1 <- c1 +c2*(d-1)
  # p4 <- dchisq(cv, df= d, ncp = delta_sq)*cv*p4.1*b
  p4.1 <- c2*(cv/2 + 1 -.75*d)+1
  p4 <-  -dchisq(cv, df= d, ncp = delta_sq)*cv*p4.1*b 
  
  p5 <- (b*big_T)^(-1)*dchisq(cv, df = d, ncp = delta_sq)*cv*g_q*w_q
  
  # Sun like: WE ARE DOING THIS ONE!! 
  if(sun_like){
    p2 <-  -delta_sq*dchisq(cv, df= (d+2), ncp = delta_sq)*cv*c2*b/2
    p4.1 <- c2*(cv/2 + 1 -.75*d)+1
    p4 <-  -dchisq(cv, df= d, ncp = delta_sq)*cv*p4.1*b 
    p4 <- 0 
  }
  
  if(b ==0){
    p5 <- 0 
  }
  
  # Type 2 error: Fail to reject null, null is false 
  the_e2 <- (p1 -p2 - p3 -p4-p5)
  
  
  return(c(the_e2))
}
