# This module contains auxiliary functions for the Application
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
create_results_table <- function(Lasso_model, #results of CV Lasso
                                 rel_Lasso_model, #results of CV relaxed Lasso
                                 model.vsurf #results of VSURF
){
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
                   "P60" = "Primary School Enrollment Rate",
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
                          "CDF" = "",
                          #"Lasso1" = NaN,
                          "Lasso" = "",
                          #"relaxed Lasso1" = NaN,
                          "relaxed Lasso" = "",
                          "Random Forest" = "")
  
  # Filling CDF
  appli_res$CDF[1:length(CDF$Var)] = CDF$unlist.CDF
  
  #Filling RF (with Rank)
  for (i in RF_vars$Var){
    appli_res$Random.Forest[appli_res$Selected.Variables == i] <- RF_vars$Rank[RF_vars$Var == i]
  }
  
  #Filling Lasso
  for (i in Lasso_coef$Var){
    #appli_res$Lasso1[appli_res$Selected.Variables == i] <- Lasso_coef$lambda.min[Lasso_coef$Var == i] #coefficient
    appli_res$Lasso[appli_res$Selected.Variables == i] <- Lasso_coef$Rank[Lasso_coef$Var == i] #rank
  }
  
  #Filling relaxed Lasso
  for (i in rel_Lasso_coef$Var){
    #appli_res$relaxed.Lasso1[appli_res$Selected.Variables == i] <- rel_Lasso_coef$lambda.min[rel_Lasso_coef$Var == i] #coefficient
    appli_res$relaxed.Lasso[appli_res$Selected.Variables == i] <- rel_Lasso_coef$Rank[rel_Lasso_coef$Var == i] #rank
  }
  
  #Creating personalised table using kable (HTML)
  options(knitr.kable.NA = "") #don't display NA
  table <- kable(appli_res,
                 col.names = gsub("[.]", " ", names(appli_res)), #proper column names
                 align = "lcccc")%>% #allignment
    kable_classic() %>%
    row_spec(c(2,4), bold = T, color = "white", background = "darkgreen") %>%
    row_spec(c(11, 12), bold = T, color = "white", background = "green") 
  
  # display results
  display_html(head(table))
}