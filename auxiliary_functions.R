# This module contains auxiliary functions for XY
# which are used in the main notebook.
# --------------------------------------------------------

# Imports
library(MASS) 
library(stats)
install.packages("glmnet")
library(glmnet) #for LASSO
install.packages("VSURF")
library(VSURF)#for RF

# Generating Sample
simulate <- function(n, #number of observations
                     p, #number of covariates
                     rho, #degree of covariance
                     beta, #vetctor of true coefficients
                     SNR # desired Signal-to-Noise ratio
){
  if (length(beta) != p){
    cat("Number of beta coefficient unequal to p")
  }else{
    #Mean of explanatory variables
    mu = rep(0,p) #all covariates are standardized with mean zero
    
    #Variance-Covariance Matrix
    ###Note: Matrix only depends on p and rho
    toep = rho^(0:(p-1)) #creates geometric series starting at one
    Sigma = toeplitz(toep) #creates toeplitz matrix from geometric series: rho^(i-j)
    
    #explanatory variables
    X = mvrnorm(n, mu, Sigma)
    
    # Set sigma based on sample variance on infinitely large test set
    var_mu = as.numeric(t(beta) %*% Sigma %*% beta)
    sigma = as.numeric(sqrt(var_mu/SNR))
    
    # Generate response variable
    Y = as.numeric(X %*% beta + rnorm(n)*sigma)
    
    #-------Creating data frame
    df <- data.frame(Y, X)
    
    list_1 = list("df" = df, "sigma" = sigma)
    return(list_1)
  } 
}

# Generate beta_1 vector
beta_1 <- function(p, #number of covariates
                   s  #degree of sparsity
){
  #--------------------------
  # Beta vector with equally 
  # spaced ones, the rest zeros
  # -------------------------
  beta = rep(0, p)
  loc = round(seq(1,p, length.out=s))
  beta[loc] = 1
  return(beta)
}

# Generate beta_2 vector
beta_2 <- function(p, #number of covariates
                   s  #degree of sparsity
){
  #--------------------------
  # Beta vector with ones at the first 
  # entries, the rest zeros
  # -------------------------
  beta = c(rep(1,s),rep(0, p-s))
  return(beta)
}

# Generate beta_3 vector
beta_3 <- function(p, #number of covariates
                   s, #degree of sparsity
                   value  #value of coefficient
){
  #--------------------------
  # Beta vector with weak sparsity:
  # first s equal to ones, the rest
  # geometric decay of values
  # -------------------------
  beta = c(rep(1,s), value^(1:(p-s))) 
  return(beta)
}

# Count retention of variables
var_retention <- function(model_coef, #coefficients of the estimated model
                          beta #true beta vector
){
  #--------------------------
  # Counts how many variables were
  # correctly identified by the estimated model
  # -------------------------
  non_zero = model_coef != 0 #find non-zero coefficients
  non_zero = as.numeric(non_zero) #transform to numeric vector
  if (length(model_coef) != length(beta)){
    non_zero = as.numeric(non_zero)[-1] # exclude intercept placeholder for lasso!
  }
  
  preserved = (non_zero + beta) == 2 # retention means double incidence
  retention = sum(preserved) # preserved is boulian vector
  return(retention)   
}

# Perform Cross-validated Lasso
cv.lasso <- function(data, #data frame - dependent variable first
                     beta # true coefficients
){
  #--------------------------
  # Uses 10 fold CV and uses 1 SE lambda
  # as conservative estimate for variable selection
  # -------------------------
  x <- data.matrix(df[,-1]) #explan var, glmnet can't use dataframe
  y <- data.matrix(df[,1]) #dependent var, glmnet can't use dataframe
  
  cv.out = cv.glmnet(x, y, alpha = 1, intercept=FALSE) # Fit lasso model on training data
  lam = cv.out$lambda.1se # Select more conservative lambda for variable selection
  
  lasso_coef = predict(cv.out, type = "coefficients", s = lam) # Display coefficients using lambda chosen by CV
  return(var_retention(lasso_coef, beta))
}



beta = beta_1(p=100,s=5)
df <- simulate(n=100, p=100, rho=0.5, beta=beta, SNR = 1)$df
x <- data.matrix(df[,-1]) #explan var, glmnet can't use dataframe
y <- data.matrix(df[,1]) #dependent var, glmnet can't use dataframe

practice.vsurf <- VSURF(x=x, y=y, mtry=100)

summary(practice.vsurf)

plot(practice.vsurf)
