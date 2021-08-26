# This module contains auxiliary functions for XY
# which are used in the main notebook.
# --------------------------------------------------------

# Imports
library(MASS) 
library(stats)
library(Matrix)
library(parallel)
library(glmnet) #for LASSO
library(VSURF)#for RF
library(ggplot2) #for plotting
library(dplyr) #for plotting
library(cowplot) #for plotting

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

# Count retention of significant variables
var_retention <- function(model_coef, #coefficients of the estimated model
                          beta #true beta vector
){
  #--------------------------
  # Counts how many significant variables were
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


# Identification of total variables
var_identification <- function(model_coef, #coefficients of the estimated model
                          beta #true beta vector
){
  #--------------------------
  # Counts how many variables were
  # correctly identified by the estimated model
  # -------------------------
  binary = model_coef != 0 #transform to binary
  binary = as.numeric(binary) #transform to numeric vector
  if (length(model_coef) != length(beta)){
    binary = as.numeric(binary)[-1] # exclude intercept placeholder for lasso!
  }
  
  corr_identified = binary == beta 
  identification = sum(corr_identified) # corr_identified is boulian vector
  return(identification)   
}

#Non-zero
var_nonzero <- function(model_coef, #coefficients of the estimated model
                               beta #true beta vector
){
  #--------------------------
  # Counts how many variables were
  # estimated as nonzero
  # -------------------------
  binary = model_coef != 0 #transform to binary
  binary = as.numeric(binary) #transform to numeric vector
  if (length(model_coef) != length(beta)){
    binary = as.numeric(binary)[-1] # exclude intercept placeholder for lasso!
  }
  
  nonzero = sum(binary)
  return(nonzero)   
}


# Perform Cross-validated Lasso
cv.lasso <- function(data, #data frame - dependent variable first
                     beta # true coefficients
){
  #--------------------------
  # Uses 10 fold CV and uses 1 SE lambda
  # as conservative estimate for variable selection
  # -------------------------
  x <- data.matrix(data[,-1]) #explan var, glmnet can't use dataframe
  y <- data.matrix(data[,1]) #dependent var, glmnet can't use dataframe
  
  cv.out = cv.glmnet(x, y, alpha = 1, intercept=FALSE) # Fit lasso model on training data
  lam = cv.out$lambda.1se # Select more conservative lambda for variable selection
  
  #---------------------
  # Retention Frequency
  #---------------------
  lasso_coef = predict(cv.out, type = "coefficients", s = lam) # Display coefficients using lambda chosen by CV
  retention = var_retention(lasso_coef, beta) #counts significant vars
  identification = var_identification(lasso_coef, beta) #counts all vars
 
  #---------------------
  # Number Nonzero elements
  #--------------------- 
  nonzero = var_nonzero(lasso_coef, beta) #count nonzero vars
  
  #---------------------
  # MSE
  #---------------------
  mse <- cv.out$cvm[cv.out$lambda == cv.out$lambda.1se] #1 standard deviation from minimum
  
  results = list("retention" = retention, "identification" =identification, "mse" = mse, "nonzero" = nonzero)
  return(results)
}

# Find retention frequency
retention_frequency <- function(results,
                                sparsity){
  res = data.frame(results)
  mean_res = as.numeric(colMeans(res)) #create list of mean values
  frequency = mean_res / true_sparsity * 100 # get percentage
  return(frequency)
}

# Find error rate
error_rate <- function(results){
  res = data.frame(results)
  loc = res == Inf # VSURF return Inf for OOB error of 0 model
  res[loc] = NA # #set INF to NA
  mean_res = as.numeric(colMeans(res, na.rm=TRUE)) #create list of mean values, ignoring NA
  return(mean_res)
}

cv.lasso_2 <- function(data, #data frame - dependent variable first
                       beta # true coefficients
){
  #--------------------------
  # Uses 10 fold CV and uses best prediciton lambda
  # as estimate for variable selection
  # -------------------------
  x <- data.matrix(data[,-1]) #explan var, glmnet can't use dataframe
  y <- data.matrix(data[,1]) #dependent var, glmnet can't use dataframe
  
  cv.out = cv.glmnet(x, y, alpha = 1, intercept=FALSE) # Fit lasso model on training data
  #lam = cv.out$lambda.1se # Select more conservative lambda for variable selection
  lam = cv.out$lambda.min
  
  #---------------------
  # Retention Frequency
  #---------------------
  lasso_coef = predict(cv.out, type = "coefficients", s = lam) # Display coefficients using lambda chosen by CV
  retention = var_retention(lasso_coef, beta) #counts significant vars
  identification = var_identification(lasso_coef, beta) #counts all vars
  
  #---------------------
  # Number Nonzero elements
  #--------------------- 
  nonzero = var_nonzero(lasso_coef, beta) #count nonzero vars
  
  #---------------------
  # MSE
  #---------------------
  mse <- cv.out$cvm[cv.out$lambda == cv.out$lambda.1se]
  
  
  results = list("retention" = retention, "identification" =identification, "mse" = mse, "nonzero" = nonzero)
  return(results)
}

#relaxed lasso
cv.relaxed_lasso <- function(data, #data frame - dependent variable first
                             beta # true coefficients
){
  #--------------------------
  # Uses 10 fold CV and uses lambda
  # and gamma minimizing prediction error
  # for variable selection
  # -------------------------
  x <- data.matrix(data[,-1]) #explan var, glmnet can't use dataframe
  y <- data.matrix(data[,1]) #dependent var, glmnet can't use dataframe
  
  cv.out = cv.glmnet(x, y,intercept=FALSE, relax=TRUE) # Fit lasso model on training data
  
  #---------------------
  # Retention Frequency
  #---------------------
  lasso_coef = predict(cv.out, type = "coefficients", s = "lambda.min", gamma = "gamma.min")#"gamma.min") # Display coefficients using lambda chosen by CV
  retention = var_retention(lasso_coef, beta) #counts significant vars
  identification = var_identification(lasso_coef, beta) #counts all vars
  
  #---------------------
  # Number Nonzero elements
  #--------------------- 
  nonzero = var_nonzero(lasso_coef, beta) #count nonzero vars
  
  #---------------------
  # MSE
  #---------------------
  mse <- cv.out$cvm[cv.out$lambda == cv.out$lambda.1se] # which gamma value is this?
  
  
  results = list("retention" = retention, "identification" =identification, "mse" = mse, "nonzero" = nonzero)
  return(results)
  
}

# RF
RF_VSURF <- function(data, #data frame - dependent variable first
                     beta #true coefficients
){
  #--------------------------
  # Uses VSURF prediction under parallelization and
  # returns number of correctly identified  significant variables.
  # Mytree and ntree are set to default
  # ------------------------- 
  x <- data.matrix(data[,-1]) #explan var, glmnet can't use dataframe
  y <- data.matrix(data[,1]) #dependent var, glmnet can't use dataframe
  
  defaultW <- getOption("warn")  #Turn off warning messages
  options(warn = -1) 
  
  
  #Variable Selection using Random Forest
  model.vsurf <- VSURF(x=x, y=y, parallel = TRUE , ncores= 4)
  
  options(warn = defaultW) #re-enable warning messages
  
  #---------------------
  # Retention Frequency
  #---------------------
  #Create boolian vector of selected coefficients
  loc = model.vsurf$varselect.pred # location of significant coefficients
  estim_var = rep(0, length(beta)) #create zero vector of correct length
  estim_var[loc] = 1 #populate zero vector
  
  retention = var_retention(estim_var, beta) #counts only significant variables
  identification = var_identification(estim_var, beta) #counts all vars
  
  #---------------------
  # Number Nonzero elements
  #--------------------- 
  nonzero = var_nonzero(estim_var, beta) #count nonzero vars
  
  #---------------------
  # OOB error
  #---------------------
  OOB_error = min(model.vsurf$err.pred) # VSURF returns a vector, final element is (minimal) OOB error and includes all !prediction! variables
  
  result = list("retention" = retention, "identification" = identification, "OOB_error" = OOB_error, "nonzero" = nonzero)
  
  return(result)
}

# Function for Violin Split
GeomSplitViolin <- ggproto(
  "GeomSplitViolin", 
  GeomViolin, 
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data, 
                      xminv = x - violinwidth * (x - xmin), 
                      xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1,'group']
    newdata <- dplyr::arrange(
      transform(data, x = if(grp%%2==1) xminv else xmaxv), 
      if(grp%%2==1) y else -y
    )
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin", 
                       grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)
# https://debruine.github.io/posts/plot-comparison/
geom_split_violin <- function (mapping = NULL, 
                               data = NULL, 
                               stat = "ydensity", 
                               position = "identity", ..., 
                               draw_quantiles = NULL, 
                               trim = TRUE, 
                               scale = "area", 
                               na.rm = FALSE, 
                               show.legend = NA, 
                               inherit.aes = TRUE) {
  layer(data = data, 
        mapping = mapping, 
        stat = stat, 
        geom = GeomSplitViolin, 
        position = position, 
        show.legend = show.legend, 
        inherit.aes = inherit.aes, 
        params = list(trim = trim, 
                      scale = scale, 
                      draw_quantiles = draw_quantiles, 
                      na.rm = na.rm, ...)
  )
}


# Function for plotting simulation results
plot_simulation_results <-function(data, #simulation data
                        beta #true beta vector
){
  #calculate true sparsity
  true_sparsity = sum(beta) #Need to change for beta type 3
  
  #Prepare the data
  ## Data frame with means
  df <- data %>% #create new data frame
    na_if(Inf) %>% #Change INF values produced by RF
    mutate(Retention = (Retention/true_sparsity)*100) %>%
    group_by(SNR, Method) %>%  
    summarize(Mean_Ret = mean(Retention, na.rm=TRUE),
              Mean_Zero = mean(Nonzero, na.rm=TRUE),
              Mean_Pred = mean(Prediction, na.rm=TRUE),
              SD_Ret = sd(Retention, na.rm= TRUE),
              SD_Zero = sd(Nonzero, na.rm=TRUE),
              SD_Pred = sd(Prediction, na.rm=TRUE))
  
  # Data frame for violin plot
  violin_data = data[data$SNR==data$SNR[600] | data$SNR==data$SNR[2000],] #selecting 0.09 and 1.22
  violin_data$SNR = as.factor(round(violin_data$SNR, digits=2))
  
  # Line breaks for SNR (logarithmic scale)
  snr.breaks = round(exp(seq(from=min(log(data$SNR)),
                             to=max(log(data$SNR)),length=4)),2)
  
  # Retention Plot
  p1 <- ggplot(data=df, aes(x=SNR, y=Mean_Ret, color=Method)) +
    geom_line(lwd=1) +
    geom_point(pch=19) +
    theme_bw() +
    xlab("Signal-to-noise ratio") +
    ylab("Retention Frequency") +
    geom_errorbar(aes(ymin=Mean_Ret-SD_Ret, ymax=Mean_Ret+SD_Ret), width=.2,
                  position=position_dodge(0.05)) +
    scale_x_continuous(trans="log", breaks=snr.breaks)
  
  # Nonzero Plot
  p2 <- ggplot(data=df, aes(x=SNR, y=Mean_Zero, color=Method)) +
    geom_line(lwd=1) +
    geom_point(pch=19) +
    theme_bw() +
    xlab("Signal-to-noise ratio") +
    ylab("Number of Nonzero Coeff") +
    geom_line(aes(x=SNR, y=true_sparsity), lwd=0.5, linetype=3, color="black") +
    geom_errorbar(aes(ymin=Mean_Zero-SD_Zero, ymax=Mean_Zero+SD_Zero), width=.2,
                  position=position_dodge(0.05)) +
    scale_x_continuous(trans="log", breaks=snr.breaks)
  
  # Plot of Prediction Error
  p3 <- ggplot(data=df, aes(x=SNR, y=Mean_Pred, color=Method)) +
    geom_line(lwd=1) +
    geom_point(pch=19) +
    theme_bw() +
    #facet_grid(rows = vars(Method)) +
    #facet_grid(formula(paste(1,"~",2))) +
    xlab("Signal-to-noise ratio") +
    ylab("Mean-squared Prediction Error") +
    geom_errorbar(aes(ymin=Mean_Pred-SD_Pred, ymax=Mean_Pred+SD_Pred), width=.2,
                  position=position_dodge(0.05)) + 
    scale_x_continuous(trans="log", breaks=snr.breaks)
  
  #Violin Plot of Nonzero distribution for selected SNR values
  p4 <- ggplot(data=violin_data, aes(x=Method, y=Nonzero, fill=SNR)) +
    geom_split_violin(color="white", trim=FALSE) +
    scale_fill_brewer(palette="Dark2") +
    theme_bw() +
    theme(legend.position=c(0.9,.75))
  
  #Start assembling plot
  #---------------------
  
  # get legend manually
  legend <- get_legend(
    # create some space to the left of the legend
    p1 + theme(legend.box.margin = margin(0, 0, 0, 12))
  )
  
  # create first grid
  g1 <- plot_grid(p1 +   theme(legend.position="none"),
                  p2 +   theme(legend.position="none"),
                  p3 +   theme(legend.position="none"),
                  p4,
                  ncol=2,
                  nrow=2)
  
  # now add the title
  title <- ggdraw() + 
    draw_label(
      "Simulation 1: n=100, p=50, Beta-type = 1, s=5, rho = 0.5",
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  
  #Create second grid with title
  g2 <- plot_grid(
    title, g1,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  )
  
  # get legend manually
  legend <- get_legend(
    # create some space to the left of the legend
    p1 + theme(legend.box.margin = margin(0, 0, 0, 12))
  )
  
  #Create final ensemble
  plot_grid(g2, legend, rel_widths = c(3, .6))
}

#-------------------
# Application
#-------------------

#Importing Salai-I-Martin 1997 data
import_millions_data <- function(){
  # Loading the data
  code_book <- read_excel("data/Salai-i-Martin_1997_data/millions.XLS", sheet=1)
  millions <- read_excel("data/Salai-i-Martin_1997_data/millions.XLS", sheet=2)
  
  # Proper column names
  colnames(code_book) = c("#", "Var_name")
  colnames(millions)[5:65] = code_book$Var_name
  
  defaultW <- getOption("warn")  #Turn off warning messages (coering empty to NA)
  options(warn = -1) 
  
  # Standardizing columns (needed for shrinkage methods)
  st_millions <- cbind(millions[,2:4], scale(millions[,5:65], center = FALSE)) #center around zero so that 0 is mean
  st_millions$gamma <- as.numeric(st_millions$gamma)
  st_millions <- st_millions %>% filter(!is.na(gamma))
  
  options(warn = defaultW) #re-enable warning messages
  
  # replace NA with mean
  colnames <- colnames(st_millions[, 3:64])
  for (i in colnames){
    mean = mean(st_millions[[i]], na.rm=TRUE)
    st_millions[[i]][is.na(st_millions[[i]])] <- mean
  }
  
  return(st_millions)
}


# Creating table of results
create_results_table <- function(Lasso_model,
                                 rel_Lasso_model,
                                 model.vsurf){
  #-------------------------------------------------------
  # This function creates a table listing
  # the important variables chosen by the three methods,
  # plus the results from Sala-I-Martin (1997)
  #------------------------------------------------------- 
  
  # 1. Retrieve Results
  #----------------------
  
  #Lasso
  Lasso_coef <- predict(Lasso_model, type="coefficients", s="lambda.min")
  Lasso_coef <- data.frame(data.matrix(Lasso_coef)) # create data frame
  Lasso_coef$Var <- row.names(Lasso_coef) # retrieve variable names
  Lasso_coef <- Lasso_coef[Lasso_coef$lambda.min != 0,]
  Lasso_coef <- Lasso_coef[sort(abs(Lasso_coef$lambda.min),decreasing=T,index.return=T)[[2]],]
  Lasso_coef$Rank <- seq(1, length(Lasso_coef$Var), by=1)
  
  
  #relaxed Lasso
  rel_Lasso_coef <- predict(rel_Lasso_model, type="coefficients", s="lambda.min")
  rel_Lasso_coef <- data.frame(data.matrix(rel_Lasso_coef)) # create data frame
  rel_Lasso_coef$Var <- row.names(rel_Lasso_coef) # retrieve variable names
  rel_Lasso_coef <- rel_Lasso_coef[rel_Lasso_coef$lambda.min != 0,]
  rel_Lasso_coef <- rel_Lasso_coef[sort(abs(rel_Lasso_coef$lambda.min),decreasing=T,index.return=T)[[2]],]
  rel_Lasso_coef$Rank <- seq(1, length(rel_Lasso_coef$Var), by=1)
  
  
  #RF
  loc = model.vsurf$varselect.pred
  RF_vars <- colnames(millions[,4:64])[loc]
  RF_vars <- data.frame(RF_vars)
  RF_vars$Rank <- row.names(RF_vars) # retrieve variable names
  colnames(RF_vars) = c("Var", "Rank")
  
  #Sala-I-Martin
  CDF <- list("GDP level 1960" = 1,
              "Fraction Confucian" = 1,
              "Life expectancy" = 0.942,
              "Equipment investment" = 0.997,
              "Sub-Saharan dummy" = 1,
              "Fraction Muslim" = 1,
              "Rule of law" = 1,
              "Open Economy" = 1,
              "Degree of capitalism"=0.987,
              "Fraction Protestant"=0.966)
  CDF <- data.frame(unlist(CDF))
  CDF$Var <- row.names(CDF) # retrieve variable names
  
  
  
  # 2. Translate column names
  #----------------------
  codebook <- list("YrsOpen" = "Open Economcy",
                   "CONFUC"="Fraction Confucian",
                   "EQINV" = "Equipment investment",
                   "LIFEE060" = "Life expectancy",
                   "GDPSH60" = "GDP level 1960",
                   "NONEQINV" = "Non-Equipment Investment",
                   "P60" = "Primary School Entrollment Rate",
                   "prightsb" = "Political rights",
                   "BUDDHA" = "Fraction Buddhist",
                   "ABSLATIT" = "Absolute Latitude",
                   "CONFUC" = "Fraction Confucian",
                   "BMS6087" = "Black Market Premium",
                   "revcoup" = "Revolution and Coups",
                   "RERD" = "Rate of Exchange Distortions",
                   "lly1" = "Ratio of Liquid Liabilities",
                   "PROT" = "Fraction Protestant")
  
  names <- names(codebook)
  for (i in names){
    Lasso_coef$Var[Lasso_coef$Var == i] = codebook[[i]]
    rel_Lasso_coef$Var[rel_Lasso_coef$Var == i] = codebook[[i]]
    RF_vars$Var[RF_vars$Var == i] = codebook[[i]]
  }
  
  # 3. Create common table
  #----------------------
  
  column1 <- unique(c(CDF$Var,Lasso_coef$Var, rel_Lasso_coef$Var, RF_vars$Var))
  appli_res <- data.frame("Selected Variables" = column1,
                          "CDF" = NaN,
                          "Lasso1" = NaN,
                          "Lasso2" = NaN,
                          "relaxed Lasso1" = NaN,
                          "relaxed Lasso2" = NaN,
                          "Random Forest" = NaN)
  
  # Filling CDF
  appli_res$CDF[1:length(CDF$Var)] = CDF$unlist.CDF
  
  #Filling RF (with Rank)
  for (i in RF_vars$Var){
    appli_res$Random.Forest[appli_res$Selected.Variables == i] <- RF_vars$Rank[RF_vars$Var == i]
  }
  
  #Filling Lasso
  for (i in Lasso_coef$Var){
    appli_res$Lasso1[appli_res$Selected.Variables == i] <- Lasso_coef$lambda.min[Lasso_coef$Var == i] #coefficient
    appli_res$Lasso2[appli_res$Selected.Variables == i] <- Lasso_coef$Rank[Lasso_coef$Var == i] #rank
  }
  
  #Filling relaxed Lasso
  for (i in rel_Lasso_coef$Var){
    appli_res$relaxed.Lasso1[appli_res$Selected.Variables == i] <- rel_Lasso_coef$lambda.min[rel_Lasso_coef$Var == i] #coefficient
    appli_res$relaxed.Lasso2[appli_res$Selected.Variables == i] <- rel_Lasso_coef$Rank[rel_Lasso_coef$Var == i] #rank
  }
  
  return(appli_res)
}