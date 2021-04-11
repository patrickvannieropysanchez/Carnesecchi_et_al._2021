# Load required packages ####
library (ggplot2)
library(robustbase)
library (minpack.lm)
library(tsoutliers)
library(investr)

# Directory setup ####
workingDir <- "/Users/patrick/Documents/PhD Local/Microscopy/FRAP2020"
setwd(workingDir)

# Read keys and Data
keyDf <- read.delim("keys.txt") # In the format of 1. file name, 2. condition

setwd(paste(workingDir, "/R_data", sep = ""))
txtList <- list.files(pattern="*.txt")

# Function calculating FRAP from acquired image intensities ####
FRAPCalc <- function(keyDf, txt){
  NucID <- gsub(".*?([0-9]+).*", "\\1", txt) # Get nuclear ID from csv file name
  
  rawData <- lapply(txt, read.delim)
  rawDf <- data.frame(rawData)
  
  blARatio <- (rawDf[("BleachArea")]/rawDf["NucleusArea"]) # Ratio area bleach relative to area whole nucleus
  FRAP <- (rawDf["BleachMean"]-rawDf["BGMean"])/(rawDf["NucleusMean"]-rawDf["BGMean"]) # Mean intensity bleach roi - mean intensity whole nucleus, both background substracted
  meanPre <- mean(FRAP[c(1,10),]) # Get average intensity pre bleaching
  FRAPNorm <- FRAP/meanPre # Performs FRAP normalisation by dividing values by average pre-bleach FRAP
  FRAPfNorm <- (FRAPNorm-FRAPNorm[11,])/(1-FRAPNorm[11,])
  FRAPfNorm2 <- (FRAPNorm-min(FRAPNorm))/(1-min(FRAPNorm)) # Performs full normalisation; 0-bleach depth, 1-recovery plateau
  
  Time2 <- rawDf["Time"]-rawDf[11, "Time"] # Sets bleaching point as t = 0
  
  resultsDf <- data.frame(c(NucID, Time2, FRAP, meanPre, FRAPNorm, FRAPfNorm,
                            FRAPfNorm2, blARatio))
  
  names(resultsDf) <- c("NucID", "Time2", "FRAP", "meanPre", "FRAPNorm", "FRAPfNorm",
                        "FRAPfNorm2", "BlARatio")
  
  calcDf <- data.frame(c(rawDf, resultsDf)) # Concatenate raw data with calculated results
  
  calcDf <- merge(calcDf, keyDf) # Merges calculated values with NucID key values
  
  return(calcDf)
}

# Model fitting functions ####
LMFitFRAP <- function(Df, modelType = 1){
  
  if (modelType == 1){
    fit <- nlrob(FRAPfNorm2 ~ A-A*exp(-Time2*k1),
                 data = Df,
                 lower = c(A = 0, k1 = 0),
                 upper = c(A = 1.5, k1 = 500),
                 start = c(A = 0.1, k1 = 0.01),
                 method = "M",
                 maxit = 1000)
  }
  
  if (modelType == 2){
    fit <- nlsLM(FRAPfNorm2 ~ (A+B)-A*exp(-Time2*k1)-B*exp(-Time2*k2), 
                 data = Df,
                 lower = c(A = 0, B = 0, k1 = 0, k2 = 0),
                 upper = c(A = 1.5, B = 1.5, k1 = 500, k2 = 500),
                 start = c(A = 0.1, B = 0.1, k1 = 0.01, k2 = 0.01)
                 )
  }
  
  return(fit)
  print(coef(fit))
}

# Fit models and export curves and coefficients #### 
PDFPlots <- function(file, Df){
  
  fitDf <- data.frame()
  IDfit <- data.frame()
  
  pdf(file, width = 11, height = 8)
  
  for (ID in unique(Df$NucID)){
    print(ID)
    
    try({
      
      model <- LMFitFRAP(Df[Df$NucID == ID &
                              Df$Time2 >= 0, ], 1)
      print(coef(model))
      
      model2 <- LMFitFRAP(Df[Df$NucID == ID &
                               Df$Time2 >= 0, ], 2)
      print(coef(model2))
      
      if (sigma(model) <= 0.1){
        
        A <- coef(model2)[["A"]]
        B <- coef(model2)[["B"]]
        
        t_half1 <- -log(0.5)/coef(model)["k1"]
        t_half2 <- invest(model2, y0 = 0.5*(A+B), x0.name = Time2, level = 0.95)[["estimate"]]
        
        AIC1 <- AIC(model, model2)["model", "AIC"]
        AIC2 <- AIC(model, model2)["model2", "AIC"]
        
        BIC1 <- BIC(model, model2)["model", "BIC"]
        BIC2 <- BIC(model, model2)["model2", "BIC"]
        
        IDfit <- cbind(ID,
                       unique(Df[Df$NucID == ID, ]$Condition),
                       coef(model)["A"],
                       coef(model)["k1"],
                       t_half1,
                       AIC1,
                       BIC1,
                       "",
                       coef(model2)["A"],
                       coef(model2)["B"],
                       coef(model2)["k1"],
                       coef(model2)["k2"],
                       t_half2,
                       AIC2,
                       BIC2,
                       "",
                       if(AIC1<AIC2) "Single" else "Double",
                       if(BIC1<BIC2) "Single" else "Double"
        )
        
        fitDf <- rbind(fitDf, IDfit)
        
        p <- ggplot(Df[Df$NucID == ID &
                         Df$Time2 >= 0, ], 
                    aes(Time2, FRAPfNorm2))+
          
          geom_point()+
          
          geom_line(aes(x = Time2, y = predict(model)), color = "blue", size = 1)+
          geom_line(aes(x = Time2, y = predict(model2)), color = "red", size = 1)+
          
          scale_x_continuous(name = "t(s)", limits = c(0, max(Df$Time2)), 
                             breaks = seq(-5, max(Df$Time2), 30), minor_breaks = NULL)+
          scale_y_continuous(name = "FRAP(t)", limits = c(0, 1.5), 
                             breaks = seq(0, 1.5, 0.1), minor_breaks = NULL)+
          
          ggtitle(paste("FRAP", ID, unique(Df[Df$NucID == ID, ]$Condition)))
        
        print(p)
      }
    })
  }
  
  colnames(fitDf) <- c("NucID",
                       "Condition",
                       "A",
                       "k1",
                       "t_half",
                       "AIC1",
                       "BIC1",
                       "",
                       "A",
                       "B",
                       "k1",
                       "k2",
                       "t_half",
                       "AIC2",
                       "BIC2",
                       "",
                       "AIC_Best",
                       "BIC_Best"
                       )
  
  write.table(fitDf, file = "/Users/patrick/Documents/PhD Local/Microscopy/FRAP2020/coefs_single_double.csv",
              sep = ",",
              row.names = F,
              col.names = T)
  
  dev.off()
}

# Create dataframe from key and data merge ####
completeDf <- data.frame()
for (txt in txtList){
  nucDf <- FRAPCalc(keyDf, txt)
    if(nrow(nucDf) >= 100){
     completeDf <- rbind(completeDf, nucDf)
  }
}

# Generate data ####
PDFPlots(("/Users/patrick/Documents/PhD Local/Microscopy/FRAP2020/R_data/FRAP_curves_single_double.pdf"), completeDf)
dev.off()
