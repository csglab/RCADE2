#  Copyright 2014 Hamed S. Najafabadi
#  
#  ********************************************************************
#  
#  This file is part of ZifRC package.
#  
#  ZifRC is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  ZifRC is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with ZifRC.  If not, see <http://www.gnu.org/licenses/>.
#  
#  ********************************************************************/

rm(list = ls())

#############################################################################

#setwd("/Users/bdogan/Documents/Research/RecognitionCode2/RCADE3/src/_R")

library(randomForest)

args <- commandArgs(trailingOnly = TRUE)

inputDir <- paste0("./tmp/",args[1]) # replace this with inputDir <- "." if testing the script outside RCADE Shell
srcDir <- "./src/_R" # replace this with srcDir <- "." if testing the script outside RCADE Shell

data <- read.csv(paste0(inputDir,"/_predict.in"),sep="\t")

#srcDir <- "/Users/bdogan/Documents/Research/RecognitionCode2/RCADE3/src/_R"
#data <- read.csv("_predict.in", sep = "\t")

nofCol <- ncol(data)

RZFInclude <- matrix(F, nrow = nrow(data), ncol = 1)
LZFInclude <- matrix(F, nrow = nrow(data), ncol = 1)
MZFInclude <- matrix(F, nrow = nrow(data), ncol = 1)

data <- cbind(data, RZFInclude, LZFInclude, MZFInclude) # add columns for logical status 

status <- matrix(0, nrow <- nrow(data), ncol = 1)

for (i in 1:nrow(data)){
  tempStr <- data[i,1]
  
  if ((regexpr('letter', tempStr)[1] == 1))
    status[i] <- 1
  else
    status[i] <- 0
  
}

diffStatus <- c(0, diff(status))

chPos1 <- which(diffStatus != 0)
chPos2 <- chPos1 - 1

chPos1 <- c(1, chPos1)
chPos2 <- c(chPos2, nrow(data))
ChangePositions <- cbind(chPos1, chPos2) # seq, letter, seq, letter, seq, letter ........

for (i in seq(from = 1, to = nrow(ChangePositions) - 1, by = 2)){
  
  Block <- data[ChangePositions[i,1]:ChangePositions[i,2],] # sequence blocks are 1, 3, 5, 7, ...., length(ChangePositions) - 1
  
  string <- Block[1,1]
  string <- as.character(string)
  indX <- which(strsplit(string, "")[[1]]=="x")
  stepInfo <- strtoi(substr(string, indX + 1, indX + 1))
  
  if (stepInfo == 2){
    
    patternRZF <- c(rep(T, stepInfo-1), F)
    patternLZF <- c(rep(F, stepInfo-1), T)
    
    seqLogicalStatusRZF <- rep(patternRZF, (nrow(Block) / stepInfo))  # logical sequence for RZF
    seqLogicalStatusLZF <- rep(patternLZF, (nrow(Block) / stepInfo))  # logical sequence for LZF
    
    data[ChangePositions[i,1]:ChangePositions[i,2], nofCol + 1] <- as.matrix(seqLogicalStatusRZF)  # add logical status of for RZF
    data[ChangePositions[i,1]:ChangePositions[i,2], nofCol + 2] <- as.matrix(seqLogicalStatusLZF)  # add logical status of for LZF
    
  }else if (stepInfo > 2){
    
    patternRZF <- c(rep(T, stepInfo-1), F)
    patternLZF <- c(F, rep(T, stepInfo-1))
    patternMZF <- c(F, rep(T, stepInfo-2), F)
    
    seqLogicalStatusRZF <- rep(patternRZF, (nrow(Block) / stepInfo))  # logical sequence for RZF
    seqLogicalStatusLZF <- rep(patternLZF, (nrow(Block) / stepInfo))  # logical sequence for LZF
    seqLogicalStatusMZF <- rep(patternMZF, (nrow(Block) / stepInfo))  # logical sequence for MZF
    
    data[ChangePositions[i,1]:ChangePositions[i,2], nofCol + 1] <- as.matrix(seqLogicalStatusRZF)  # add logical status of for RZF
    data[ChangePositions[i,1]:ChangePositions[i,2], nofCol + 2] <- as.matrix(seqLogicalStatusLZF)  # add logical status of for LZF
    data[ChangePositions[i,1]:ChangePositions[i,2], nofCol + 3] <- as.matrix(seqLogicalStatusMZF)  # add logical status of for LZF
  }
}

############################################################################################

# Biochemical properties of amino acids

ExtractFeature <- function(PSeq){
  
  biochem <- read.table( paste0(srcDir,"/biochem2d.txt"), header = F, sep = "\t")
  biochem <- as.matrix(biochem)
  
  feat <- NULL
  for (j in 1:length(PSeq)){
    aa <- PSeq[j]
    ind <- which(aa == biochem[,1])
    feat <- cbind(feat, biochem[ind, 2:ncol(biochem)])
  }
  
  featureVector <- matrix(feat, nrow = 1, ncol = (length(PSeq)*2))
  
  features <- matrix(0, nrow = 1, ncol = length(featureVector))
  for (j in 1:length(featureVector))
    features[j] <- as.numeric(featureVector[j])
  
  return(features)
}

############################################################################################

# Function for prediction
PredictbyRF <- function(InputData, model){
  
  # RF model requires column names.
  cname <- NULL
  for (i in 1:length(InputData)){
    s <- c("C", as.character(i))
    cname <- cbind(cname, paste(s, sep="", collapse="") )
  }
  
  cname <- make.names(cname)
  
  test <- InputData
  colnames(test) <- cname[1:ncol(InputData)]
  
  test <- as.data.frame(test)
  
  pred <- predict(model, test)
  
  return(pred)
  
}

############################################################################################
# Best ZF models based on number of features for each output position
# "1" use 4 residues
# "2" use 7 residues
# "3" use 12 residues
BestSZF <- c(3, 1, 3, 1, 1, 1, 2, 3, 3, 2, 3, 3)
BestMZF <- c(3, 3, 1, 3, 2, 1, 1, 2, 1, 1, 1, 3)
BestLZF <- c(3, 3, 1, 3, 3, 1, 2, 3, 2, 1, 2, 3)
BestRZF <- c(3, 3, 3, 3, 1, 1, 1, 1, 3, 2, 1, 3)

RZFcase <- c("RZF", "RZF", "RZF", "RZF", "RZF", "RZF", "RZF", "SZF", "RZF", "RZF", "RZF", "SZF") # Model choice for RZF case
LZFcase <- c("LZF", "LZF", "LZF", "LZF", "LZF", "LZF", "SZF", "SZF", "SZF", "LZF", "LZF", "SZF") # Model choice for LZF case
MZFcase <- c("LZF", "LZF", "LZF", "RZF", "RZF", "LZF", "RZF", "SZF", "RZF", "RZF", "RZF", "SZF") # Model choice for MZF case


# Read RF Models
SZFModels <- list()
RZFModels <- list()
LZFModels <- list()
MZFModels <- list()

for (i in 1:12){
  SZFModels[[i]] <- readRDS(paste0(srcDir,"/Models/SZF/", i, ".rds"))
  RZFModels[[i]] <- readRDS(paste0(srcDir,"/Models/RZF/", i, ".rds")) 
  LZFModels[[i]] <- readRDS(paste0(srcDir,"/Models/LZF/", i, ".rds")) 
  MZFModels[[i]] <- readRDS(paste0(srcDir,"/Models/MZF/", i, ".rds")) 
}
############### construct the matrix of predictors

indx <- which(status == 0)    # select real ZFs
covariates <- data[indx,2:13] # to save computational time do predictions for only real ZFs
LogStatus <- data[indx, (nofCol+1):(nofCol+3)] # ZF logical status

PredictedOutputs <- NULL

for (k in 1:nrow(covariates)){
  
  logStat <- LogStatus[k,]
  
  if (logStat[1] == FALSE && logStat[2] == FALSE && logStat[3] == FALSE){   # must use SZF table
    
    result <- NULL
    
    for (i in 1:12){ # for each output position
      
      featExtMethod <- BestSZF[i] # which features, 4 residues, 7 residues or 12 residues
      
      if (featExtMethod == 1)
        AASeq <- as.matrix(covariates[k, c(6, 8, 9, 12)])
      else if (featExtMethod == 2)
        AASeq <- as.matrix(covariates[k, c(3, 5, 6, 7, 8, 9, 12)])
      else if (featExtMethod == 3)
        AASeq <- as.matrix(covariates[k,])
      
      InputData <- ExtractFeature(AASeq)
      pred <- PredictbyRF(InputData, SZFModels[[i]])
      
      result <- cbind(result, pred)
      
    } # end for (i in 1:12)
    
    
  }else if (logStat[1] == TRUE && logStat[2] == FALSE && logStat[3] == FALSE){  # use table for RZF case
    
    result <- NULL
    
    for (i in 1:12){ # for each output position
      
      if (RZFcase[i] == "SZF"){
        
        featExtMethod <- BestSZF[i] # which features, 4 residues, 7 residues or 12 residues
        
        if (featExtMethod == 1)
          AASeq <- as.matrix(covariates[k, c(6, 8, 9, 12)])
        else if (featExtMethod == 2)
          AASeq <- as.matrix(covariates[k, c(3, 5, 6, 7, 8, 9, 12)])
        else if (featExtMethod == 3)
          AASeq <- as.matrix(covariates[k,])
        
        InputData <- ExtractFeature(AASeq)
        pred <- PredictbyRF(InputData, SZFModels[[i]])
        
        result <- cbind(result, pred)
        
      }else{
        
        featExtMethod <- BestRZF[i] # which features, 4 residues, 7 residues or 12 residues
        
        if (featExtMethod == 1)
          AASeq <- c(as.matrix(covariates[k, c(6, 8, 9, 12)]), as.matrix(covariates[k+1, c(6, 8, 9, 12)]))
        else if (featExtMethod == 2)
          AASeq <- c(as.matrix(covariates[k, c(3, 5, 6, 7, 8, 9, 12)]), as.matrix(covariates[k+1, c(3, 5, 6, 7, 8, 9, 12)]))
        else if (featExtMethod == 3)
          AASeq <- c(as.matrix(covariates[k,]), as.matrix(covariates[k+1,]))
        
        InputData <- ExtractFeature(AASeq)
        pred <- PredictbyRF(InputData, RZFModels[[i]])
        
        result <- cbind(result, pred)
      }
      
    } # end for (i in 1:12)
    
  }else if (logStat[1] == FALSE && logStat[2] == TRUE && logStat[3] == FALSE){  # use table for LZF case
    
    result <- NULL
    
    for (i in 1:12){ # for each output position
      
      if (LZFcase[i] == "SZF"){
        
        featExtMethod <- BestSZF[i] # which features, 4 residues, 7 residues or 12 residues
        
        if (featExtMethod == 1)
          AASeq <- as.matrix(covariates[k, c(6, 8, 9, 12)])
        else if (featExtMethod == 2)
          AASeq <- as.matrix(covariates[k, c(3, 5, 6, 7, 8, 9, 12)])
        else if (featExtMethod == 3)
          AASeq <- as.matrix(covariates[k,])
        
        InputData <- ExtractFeature(AASeq)
        pred <- PredictbyRF(InputData, SZFModels[[i]])
        
        result <- cbind(result, pred)
        
      }else{
        
        featExtMethod <- BestLZF[i] # which features, 4 residues, 7 residues or 12 residues
        
        if (featExtMethod == 1)
          AASeq <- c(as.matrix(covariates[k, c(6, 8, 9, 12)]), as.matrix(covariates[k-1, c(6, 8, 9, 12)]))
        else if (featExtMethod == 2)
          AASeq <- c(as.matrix(covariates[k, c(3, 5, 6, 7, 8, 9, 12)]), as.matrix(covariates[k-1, c(3, 5, 6, 7, 8, 9, 12)]))
        else if (featExtMethod == 3)
          AASeq <- c(as.matrix(covariates[k,]), as.matrix(covariates[k-1,]))
        
        InputData <- ExtractFeature(AASeq)
        pred <- PredictbyRF(InputData, LZFModels[[i]])
        
        result <- cbind(result, pred)
      }
      
    } # end for (i in 1:12)
    
  }else if (logStat[1] == TRUE && logStat[2] == TRUE && logStat[3] == TRUE){ # use table for MZF case
    
    result <- NULL
    
    for (i in 1:12){ # for each output position
      
      if (MZFcase[i] == "SZF"){
        
        featExtMethod <- BestSZF[i] # which features, 4 residues, 7 residues or 12 residues
        
        if (featExtMethod == 1)
          AASeq <- as.matrix(covariates[k, c(6, 8, 9, 12)])
        else if (featExtMethod == 2)
          AASeq <- as.matrix(covariates[k, c(3, 5, 6, 7, 8, 9, 12)])
        else if (featExtMethod == 3)
          AASeq <- as.matrix(covariates[k,])
        
        InputData <- ExtractFeature(AASeq)
        pred <- PredictbyRF(InputData, SZFModels[[i]])
        
        result <- cbind(result, pred)
        
      }else if (MZFcase[i] == "RZF"){
        
        featExtMethod <- BestRZF[i] # which features, 4 residues, 7 residues or 12 residues
        
        if (featExtMethod == 1)
          AASeq <- c(as.matrix(covariates[k, c(6, 8, 9, 12)]), as.matrix(covariates[k+1, c(6, 8, 9, 12)]))
        else if (featExtMethod == 2)
          AASeq <- c(as.matrix(covariates[k, c(3, 5, 6, 7, 8, 9, 12)]), as.matrix(covariates[k+1, c(3, 5, 6, 7, 8, 9, 12)]))
        else if (featExtMethod == 3)
          AASeq <- c(as.matrix(covariates[k,]), as.matrix(covariates[k+1,]))
        
        InputData <- ExtractFeature(AASeq)
        pred <- PredictbyRF(InputData, RZFModels[[i]])
        
        result <- cbind(result, pred)
        
        
      }else if (MZFcase[i] == "LZF"){
        
        featExtMethod <- BestLZF[i] # which features, 4 residues, 7 residues or 12 residues
        
        if (featExtMethod == 1)
          AASeq <- c(as.matrix(covariates[k, c(6, 8, 9, 12)]), as.matrix(covariates[k-1, c(6, 8, 9, 12)]))
        else if (featExtMethod == 2)
          AASeq <- c(as.matrix(covariates[k, c(3, 5, 6, 7, 8, 9, 12)]), as.matrix(covariates[k-1, c(3, 5, 6, 7, 8, 9, 12)]))
        else if (featExtMethod == 3)
          AASeq <- c(as.matrix(covariates[k,]), as.matrix(covariates[k-1,]))
        
        InputData <- ExtractFeature(AASeq)
        pred <- PredictbyRF(InputData, LZFModels[[i]])
        
        result <- cbind(result, pred)
        
        
      }else if (MZFcase[i] == "MZF"){
        
        featExtMethod <- BestMZF[i] # which features, 4 residues, 7 residues or 12 residues
        
        if (featExtMethod == 1)
          AASeq <- c(as.matrix(covariates[k-1, c(6, 8, 9, 12)]), as.matrix(covariates[k, c(6, 8, 9, 12)]), as.matrix(covariates[k+1, c(6, 8, 9, 12)]))
        else if (featExtMethod == 2)
          AASeq <- c(as.matrix(covariates[k-1, c(3, 5, 6, 7, 8, 9, 12)]), as.matrix(covariates[k, c(3, 5, 6, 7, 8, 9, 12)]), as.matrix(covariates[k+1, c(3, 5, 6, 7, 8, 9, 12)]))
        else if (featExtMethod == 3)
          AASeq <- c(as.matrix(covariates[k-1,]), as.matrix(covariates[k,]), as.matrix(covariates[k+1,]))
        
        InputData <- ExtractFeature(AASeq)
        pred <- PredictbyRF(InputData, MZFModels[[i]])
        
        result <- cbind(result, pred)
        
      }
      
    } # end for (i in 1:12)
    
  }  
  
  
  # normalize the predicted output
  prOutput <- matrix(result, nrow = 3, ncol = 4, byrow = TRUE)
  for (j in 1:3)
    prOutput[j,] <- prOutput[j,] / sum(prOutput[j,])
  
  result <- c(prOutput[1,], prOutput[2,], prOutput[3,])
  
  PredictedOutputs <- rbind(PredictedOutputs, result)  


} #end for (k in 1:nrow(covariates))


for (i in 1:12 )
  data[indx,i+13] <- PredictedOutputs[,i]


#write.table(data, file = "_predict.RF.out",sep="\t",row.names=F)


write.table(data,file=paste0(inputDir,"/_predict.RF.out"),sep="\t",row.names=F)


rm(covariates)
rm(data)


