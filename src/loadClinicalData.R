setwd("/imppc/labs/mplab/share/metallopeptidases/")

cancers <- c("CRC", "PAAD",  "KIRP", "LIHC", "BLCA", "PRAD", "LUSC", "THCA", "KIRC", "UCEC", "LUAD", "BRCA")
cancers <- sort(cancers)


for (cancer in cancers){
  
  #cancer <- "CRC"
  print(cancer)

  donordata <- read.table(paste0("icgc/",cancer,"/donor/donor.",cancer,"-US.tsv"),header = T, sep = "\t")
  rownames(donordata) <- donordata$submitted_donor_id
  
  donordata2 <- read.csv(paste0("tcga.clinical/TCGA-",cancer,"_clinical.csv"))
  rownames(donordata2) <- donordata2$submitter_id
  
  donors <- merge(donordata,donordata2,by=0,all=T)
  rownames(donors) <- donors$Row.names
  donors <- donors[,-1]
  save(donors, file=paste0("data/clinical_",cancer,".rda"))
  
}
