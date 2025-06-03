
# -------------------------------------------------------------------------
# Fitted CV ---------------------------------------------------------------
# -------------------------------------------------------------------------
fitted_model <- function(d=1, cv_matrix, alpha){
  try_b <- cv_matrix[,1]
  alpha_levels <- 1-alpha 
  chisq_cv <-  qchisq(alpha_levels, df = d)/d 
  specific_cvs <- cv_matrix[,d+1]
  
  # Get a CV fitted as a function of b 
  fit <-  lm(specific_cvs ~  poly(try_b, 3, raw = T)) 
  fit <- summary(fit)
  
  return_me <- c(fit$coefficients[,1], fit$adj.r.squared, d)
  names(return_me) <- c("intercept", "beta1", "beta2", "beta3", "adj.r.sq", "d")
  
  return(return_me)
}


# -------------------------------------------------------------------------
# Generate all Fitted CV --------------------------------------------------
# -------------------------------------------------------------------------

files <- list.files()
files <- grep("\\.csv",files, value = T)
files <- files[-which(files == "fitted_CV.csv")]

# Read in the file, and create a fitted linear regression line 
# Keep information about the file: mother kernel, lugsail, alpha, d

the_fits <- data.frame(NA, NA, NA,NA, NA, NA, NA, NA, NA)
colnames(the_fits) <- c("intercept", "beta1", "beta2", "beta3", "adj.r.sq","d",
                        "kernel", "lugsail", "alpha")
for(the_file in files){
  cv_matrix <- read.csv(the_file)
  the_file <- strsplit(the_file, "_")[[1]]
  mother_kernel <- the_file[1]
  lugsail_type  <- the_file[2]
  alpha <- as.numeric(paste(".", gsub("\\.csv","",the_file[3]), sep=""))
  fit <- fitted_model(1, cv_matrix, alpha)
  fit <- c(fit, kernel = mother_kernel , lugsail = lugsail_type, alpha = alpha)
  the_fits <- rbind(the_fits, fit)
}

# Store the results for easy access 
the_fits <- the_fits[-1, ]
the_fits$lugsail <- factor(the_fits$lugsail, 
                           levels = c("Mother", "Zero", "Adapt", "Over"), 
                           order = T)
the_fits <- the_fits[order(the_fits$kernel, the_fits$lugsail, the_fits$alpha), ]
write.csv(the_fits, "fitted_CV.csv",row.names = F) 


# -------------------------------------------------------------------------
# Find Specific CV --------------------------------------------------------
# -------------------------------------------------------------------------

# New_b is a vector of b values (could have length >1) of which 
# you want to find the CV for. 

get_cv <- function(d, alpha, the_kernel, lugsail, new_b){
  # Read in all fitted values 
  the_fits <- read.csv(fitted_cv)
  chisq_cv <-  qchisq(1-alpha, df = d)/d 
  
  # Pull out only the values you need 
  index <- the_fits$kernel == the_kernel & the_fits$lugsail == lugsail &
    the_fits$alpha == alpha & the_fits$d == d
  
  coefficients <- the_fits[index, c("beta1", "beta2","beta3")]
  intercept <-  the_fits[index, c("intercept")]
  
  # Fitted value of b 
  new_b<-data.frame(poly(new_b, 3, raw = T))
  cv_by_b <- apply(new_b, 1, function(x) sum(x*coefficients))
  cv_by_b <- intercept + cv_by_b
  
  return(cv_by_b)
}