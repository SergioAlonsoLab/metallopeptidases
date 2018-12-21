# Calculate the genomic distance between all selected regions of interest (ROI)
# by Sergio Alonso

library(GenomicRanges)
# library(bigmemory)
setwd("/imppc/labs/mplab/share/metallopeptidases/")

load("data/regions.rda",verbose=T)

# simplified version: calculate distances by chromosome

chrs <- paste0("chr",c(1:22,"X","Y"))

distances <- lapply(chrs,function(chr) {
  cat(sprintf("Calculating distances for %s\n",chr))
  roisInCHR <- ROI[seqnames(ROI) == chr]
  x <- matrix(NA,nrow=length(roisInCHR),ncol=length(roisInCHR)) # to allocate memory space
  colnames(x) <- rownames(x) <- names(roisInCHR)
  for(i in 1:ncol(x)) {
    x[i,] <- distance(roisInCHR[i],roisInCHR,ignore.strand=T) 
  }
  return(x)
})

names(distances) <- chrs

saveRDS(distances,"data/distances.rds")
