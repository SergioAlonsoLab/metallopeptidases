#Average methylation over CGI
#Sergio Alonso

library(GenomicRanges)
setwd("/imppc/labs/mplab/share/metallopeptidases/")
load("data/regions.rda", v=T)

outputdir <- "/imppc/labs/mplab/share/ICGC"

cancers <- commandArgs(trailingOnly = T)

# mean by ROI

rois <- unique(geneToProbe[,c("ROI","probeID")])

meanByROI <- function(betas) {
  shared <- rois$probeID %in% rownames(betas)
  x <- betas[rois$probeID[shared],]
  apply(x,2,function(sample) {
    tapply(sample,rois$ROI[shared],mean,na.rm=T)
  })
}


for (cancer in cancers){
  cat(cancer,"\n")
  betas <- readRDS(sprintf("%s/%s/methICGC/correctedBetas.rds",outputdir,cancer))

  cat("Calcualting means by ROI\n")
  betasROI <- meanByROI(betas)
  saveRDS(betasROI, file=sprintf("%s/%s/methICGC/betasROI.rds",outputdir,cancer))
  
  cat("Calcualting somatic changes by ROI\n")
  normals <- substr(colnames(betasROI),14,15) == "11"
  somaticROI <- betasROI - rowMeans(betasROI[,normals],na.rm=T)
  
  saveRDS(somaticROI, file=sprintf("%s/%s/methICGC/somaticROI.rds",outputdir,cancer))
  
}




