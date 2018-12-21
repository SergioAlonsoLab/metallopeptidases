# correct beta values type II
# remove samples of array HM27
# and probes not in array HM450
# by Sergio Alonso

# This script works from Rstudio, but not from Rscript (!?)

library(tools)

cancers <- c("BLCA",
             "BRCA",
             
             "CRC",
             
             "KIRC",
             "KIRP",
             "LAML",
             "LIHC",
             "LUAD",
             "LUSC",
             "PAAD",
             "PRAD",
             
             "SKCM",
             "STAD",
             "THCA",
             "UCEC")


if(length(commandArgs(trailingOnly = T)) > 0) cancers <- commandArgs(trailingOnly = T)

setwd("/imppc/labs/mplab/share/metallopeptidases/")

load("data/probeInfo.rda",v=T)
probeInfo <- as.data.frame(probeInfo)

source("/imppc/labs/mplab/share/Illumina450/src/correct_Beta_values.R")

for (cancer in cancers) {
  cat(cancer,"\n")
  setwd(paste0("/imppc/labs/mplab/share/ICGC/",cancer,"/methICGC"))
  
  load("betas.rda",v=T)
  
  cat(sprintf("%s methylation data loaded. %i rows x %i columns\n",cancer,nrow(betas),ncol(betas)))
  
  hm450 <- which(apply(is.na(betas),2,sum) < 1e5)
  validProbes <- rownames(betas)[which(apply(!is.na(betas[,hm450]),1,sum) >= 10)]
  validProbes <- validProbes[substr(validProbes,1,2) == "cg"]
  
  cat(head(validProbes),sep="\n")
  cat("...\n")
  cat(tail(validProbes),sep="\n")
  
  cat("Selecting samples and probes\n")
  betasCorrected <- betas[validProbes,hm450] # only HM450 samples
  
  cat(sprintf("Selected %i probes and %i samples\n",nrow(betasCorrected),ncol(betasCorrected)))
  
  probeType <- probeInfo[validProbes,"designType"]

  cat("Correcting betas\n")
  for(i in 1:ncol(betasCorrected)) {
    betasCorrected[,i] <- correct.Beta(betasCorrected[,i],probeType)
  }
  
  ## Combine duplicates
  
  samples <- colnames(betasCorrected) <- substr(colnames(betasCorrected),1,16)
  dups <- duplicated(samples)
  
  cat("Removing duplicates\n")
  for(sample in unique(samples[dups])) {
    selcols <- which(samples==sample)
    betasCorrected[,selcols] <- rowMeans(betasCorrected[,selcols],na.rm = T)
  }
  
  betasCorrected <- betasCorrected[,samples[!dups]]
  
  cat("Saving corrected betas\n")
  saveRDS(betasCorrected,file="correctedBetas.rds")
}  

