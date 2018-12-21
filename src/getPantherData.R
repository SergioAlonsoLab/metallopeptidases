#Get Panther DB data
setwd("/imppc/labs/mplab/share/metallopeptidases/")
library(data.table)
library("biomaRt")

ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org") #uses human ensembl annotations

attributes <- listAttributes(ensembl)
subset(attributes,grepl("hgnc",attributes$name))

pantherBM <- getBM(attributes=c('gene_biotype',
                                'hgnc_symbol',
                                'external_gene_name',
                                'hmmpanther'),
                   filters='transcript_biotype',
                   values='protein_coding',
                   mart = ensembl)

pantherBM <- pantherBM[pantherBM$hmmpanther!="",]
# pantherBM$hmmpanther <- gsub(":.*","",pantherBM$hmmpanther) # this line removed the subfamilies... not a good strategy
pantherBM <- unique(pantherBM)

#ftp://ftp.pantherdb.org/hmm_classifications/current_release
pantherdb <- fread("data/PANTHER13.1_HMM_classifications")[,1:2]
colnames(pantherdb)<-c("id","name")
# pantherdb <- pantherdb[!grepl(":", pantherdb$id),]

panther <- merge(pantherBM,pantherdb,by.x="hmmpanther", by.y="id", all.x = T,sort=F)

x1 <- data.frame(name=tapply(panther$name,panther$hmmpanther,unique))
x2 <- data.frame(genes=tapply(panther$hgnc_symbol,panther$hmmpanther,function(x) {
  x %>% unique %>% sort -> x
  x[x!=""] -> x
  paste(x,collapse=";")
}))

all(rownames(x1)==rownames(x2))

panther2 <- cbind(x1,x2)

rm(x1,x2)

panther2$hmmpanther <- rownames(panther2)
panther2$Freq <- sapply(panther2$genes,function(x) length(unlist(strsplit(x,split=";"))))

familyToGene <- function(familyNames) {
  panther2$genes[panther2$hmmpanther %in% familyNames] %>% strsplit(.,split = ";") %>% unlist %>% unique
}

geneToFamily <- function(geneNames) {
  panther$hmmpanther[panther$hgnc_symbol %in% geneNames] %>% unique 
}


save(panther,panther2,familyToGene,geneToFamily,file="data/pantherdb.rda")
