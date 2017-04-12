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

#setwd("/Users/bdogan/Documents/Research/tbtkProject/nature21683/RFiles/rcadeHamed/RCADER-master/src/NewModel")

library(randomForest)

args <- commandArgs(trailingOnly = TRUE)

inputDir <- paste0("./tmp/",args[1]) # replace this with inputDir <- "." if testing the script outside RCADE Shell
srcDir <- "./src/_R" # replace this with srcDir <- "." if testing the script outside RCADE Shell

data <- read.csv(paste0(inputDir,"/_predict.in"),sep="\t")
#data <- read.csv("_predict.in", sep = "\t")

NZFInclude <- matrix(0, nrow = nrow(data), ncol = 1)

data <- cbind(data, NZFInclude) # add the logical status to decide wheter to include the neighboring ZF or not

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
  
  pattern <- c(rep("T", stepInfo-1), "F")
  
  seqLogicalStatus <- rep(pattern, (nrow(Block) / stepInfo))  # if NZF will be added or not
  
  data[ChangePositions[i,1]:ChangePositions[i,2], ncol(data)] <- as.matrix(seqLogicalStatus)  # add logical status of each ZF to data
  
}

############################################################################################

# Method 3 biochemical properties (6 property) of amino acids

ExtractFeature <- function(PSeq){
  
  biochem <- read.table( paste0(srcDir,"/biochemprop.txt"), header = F, sep = "\t")
  biochem <- as.matrix(biochem)
  
  feat <- NULL
  for (j in 1:length(PSeq)){
    aa <- PSeq[j]
    ind <- which(aa == biochem[,1])
    feat <- cbind(feat, biochem[ind, 2:ncol(biochem)])
  }

  featureVector <- matrix(feat, nrow = 1, ncol = (length(PSeq)*6))
  
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
  for (i in 1:ncol(InputData)){
    s <- c("C", as.character(i))
    cname <- cbind(cname, paste(s, sep="", collapse="") )
  }
  
  cname <- make.names(cname)
  
  Test <- NULL
  Test <- InputData
  colnames(Test) <- cname[1:ncol(InputData)]
  
  TestDF <- as.data.frame(Test)
  
  # Use one of the random forest model for each output variable.
  result <- NULL
  
  for (j in 1:12){
    
    pred <- NULL
    for ( i in 1:nrow(TestDF)){
      test = TestDF[i,]
      pred <- rbind(pred, predict(model[[j]], test)) # Prediction for the input test data
    }
    
    result <- cbind(result, pred)
  }
  
  # normalize the predicted output
  prOutput <- matrix(result, nrow = 3, ncol = 4, byrow = TRUE)
  for (j in 1:3)
    prOutput[j,] <- prOutput[j,] / sum(prOutput[j,])
  
  result <- c(prOutput[1,], prOutput[2,], prOutput[3,])

  return(result)
  
}

############################################################################################
# read models
modelNZF <- readRDS(paste0(srcDir,"/modelNZF.rds"))       # model for prediction by using neighboring ZF
modelSingle <- readRDS(paste0(srcDir,"/modelSingle.rds")) # model for regular prediction

############### construct the matrix of predictors
covariates <- data[,c(7,9,10,13)]
NZFStatus <- data[, ncol(data)] # ZF logical status

PredictedOutputs <- NULL

for (k in 1:nrow(covariates)){
  
  logStat <- NZFStatus[k]
  
  if (logStat == "T"){
    
    AAseq <- c(as.matrix(covariates[k,]), as.matrix(covariates[k+1,]))
    InputData <- ExtractFeature(AAseq) # -----------------------------------> Feature extraction by method 3
    
    result <- PredictbyRF(InputData, modelNZF)

    
  }else if (logStat == "F" || logStat == 0){
    
      AAseq <- as.matrix(covariates[k,])
      InputData <- ExtractFeature(AAseq) # -----------------------------------> Feature extraction by method 3

      result <- PredictbyRF(InputData, modelSingle)
      
  }
  
  PredictedOutputs <- rbind(PredictedOutputs, result)  
}


for (i in 1:12 )
  data[,i+13] <- PredictedOutputs[,i]

#write.table(data, file = "_predict.RF.out",sep="\t",row.names=F)

write.table(data,file=paste0(inputDir,"/_predict.RF.out"),sep="\t",row.names=F)


rm(covariates)
rm(data)


