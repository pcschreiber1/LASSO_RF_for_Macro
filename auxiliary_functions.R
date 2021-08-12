# This module contains auxiliary functions for XY
# which are used in the main notebook.
# --------------------------------------------------------

# Imports
library(MASS) 
library(stats)

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
    
    # Set snr based on sample variance on infinitely large test set
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
