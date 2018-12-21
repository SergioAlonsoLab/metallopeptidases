# get gene info from probes and define regions of interest
# Andreu Alibés & Sergio Alonso

# libraries and constants ------

library(GenomicRanges)
library(dplyr) # for readability, some consecutive functions are piped

PROMOTERWIDTH <- 2000


setwd("/imppc/labs/mplab/share/metallopeptidases/")

# Generate the data using the most updated information in the ENSEMBL: ENST GRanges object ----

library(biomaRt)

ensmblMart <- useMart("ENSEMBL_MART_ENSEMBL","hsapiens_gene_ensembl")
martAtt <- listAttributes(ensmblMart)

attrs <- c(
  "hgnc_symbol",
  "ensembl_transcript_id",
  "transcript_biotype",
  "transcript_appris",
  "transcript_tsl",
  "chromosome_name",
  "transcription_start_site",
  "transcript_start",
  "transcript_end",
  "transcript_length",
  "strand"
)

human_coding_transcripts <- getBM(attributes = attrs,
                               filters = c("transcript_biotype","chromosome_name"),
                               values = list("protein_coding",c(1:22,"X","Y","M")),
                               mart = ensmblMart) # this might take ~ 1m



ENST <- GRanges(seqnames = paste0("chr",human_coding_transcripts$chromosome_name),
        ranges=IRanges(start = human_coding_transcripts$transcript_start,
                        end = human_coding_transcripts$transcript_end),
        strand = human_coding_transcripts$strand)

names(ENST) <- human_coding_transcripts$ensembl_transcript_id
mcols(ENST) <- human_coding_transcripts

TSS <- GRanges(seqnames = paste0("chr",human_coding_transcripts$chromosome_name),
                ranges=IRanges(start = human_coding_transcripts$transcription_start_site - PROMOTERWIDTH/2,
                               end = human_coding_transcripts$transcription_start_site + PROMOTERWIDTH/2),
                strand = "*")

names(TSS) <- human_coding_transcripts$ensembl_transcript_id

# Probes as GenomicRanges: probeInfo GRanges object ----

load("data/probeInfo.rda",verbose=T)
probeInfo <- probeInfo[mcols(probeInfo)$probeType=="cg"] # only methylation probes

table(strand(probeInfo)) # 30 probes with no strand info
table(width(probeInfo),strand(probeInfo)) # and no position information

probeInfo <- subset(probeInfo,strand!="*") # remove those 30 probes
table(seqnames(probeInfo))

# CGIs & Shores as Genomic Ranges: CGI GRanges object ----

inCGI <- probeInfo$CGIposition %in% c("N_Shore","S_Shore","Island")

probeInfo$CGI[!inCGI] <- NA

CGI <- as.data.frame(probeInfo[inCGI])

coor <- tapply(CGI$start,CGI$CGI,function(x) range(x,na.rm=T))
chr <- tapply(as.character(CGI$seqnames),CGI$CGI,unique)

CGI <- GRanges(seqnames = chr,
        IRanges(start = sapply(coor,function(x) x[1] - 1),
                end = sapply(coor,function(x) x[2])) + 1)

names(CGI) <- sprintf("CGI:%s:%i-%i",seqnames(CGI),start(CGI),end(CGI))

# intersect TSS with probes, to name ROIs: pbTSS GRanges object ----

hits <- findOverlaps(TSS,probeInfo,ignore.strand=T)

tssToProbe <- data.frame(ensembl_transcript_id = names(ENST)[hits@from],
                         tss = ENST@elementMetadata$transcription_start_site[hits@from],
                         probeID = names(probeInfo)[hits@to],
                         chr = seqnames(probeInfo)[hits@to],
                         probeCpG = start(probeInfo)[hits@to],
                         stringsAsFactors = F)

x <- t(sapply(tapply(tssToProbe$probeCpG,tssToProbe$ensembl_transcript_id,range),c))

tssToProbe$TSSstart <- x[tssToProbe$ensembl_transcript_id,1]-1
tssToProbe$TSSend <- x[tssToProbe$ensembl_transcript_id,2]+1
tssToProbe$TSSlength <- tssToProbe$TSSend-tssToProbe$TSSstart

tssToProbe$TSS <- sprintf("TSS:%s:%i-%i",tssToProbe$chr,tssToProbe$TSSstart,tssToProbe$TSSend)

duplicates <- duplicated(paste(tssToProbe$ensembl_transcript_id,tssToProbe$TSS))
enstToTSS <- tssToProbe[!duplicates,c("ensembl_transcript_id","TSS","chr","TSSstart","TSSend")]


pbTSS <- GRanges(seqnames=enstToTSS$chr,
               IRanges(start=enstToTSS$TSSstart,
                       end=enstToTSS$TSSend))

names(pbTSS) <- enstToTSS$TSS

pbTSS <- unique(pbTSS) # several ENST can hit the same pbTSS

# Start creating the gene-transcript-ROI-probe-conversion data.frame ----

geneToTranscript <- ENST@elementMetadata[,c("hgnc_symbol","ensembl_transcript_id","transcription_start_site")]

# Find overlaps between TSS and CGI:  ----

hits <- findOverlaps(TSS,CGI,ignore.strand=T)

enstToCGI <- data.frame(ensembl_transcript_id = names(TSS)[hits@from],
                       CGI = names(CGI)[hits@to],
                       stringsAsFactors = F)

transcriptToROI <- merge(enstToCGI,enstToTSS[,1:2],by="ensembl_transcript_id",all.x=T,all.y=T)
transcriptToROI$ROI <- ifelse(is.na(transcriptToROI$CGI),transcriptToROI$TSS,transcriptToROI$CGI)


geneToROI <- merge(geneToTranscript,
             transcriptToROI[,c("ensembl_transcript_id","ROI")],
             by="ensembl_transcript_id",
             sort=F,
             all.x=T)

geneToROI <- as.data.frame(geneToROI)

# 10971 rows do not have an associated ROI (ENS without probes)

geneToROI <- na.omit(geneToROI)


ROI <- c(CGI,pbTSS)[unique(geneToROI$ROI)]

hits <- findOverlaps(ROI,probeInfo,ignore.strand=T)

ROItoProbe <- data.frame(ROI = names(ROI)[hits@from],
                         probeID = names(probeInfo)[hits@to],
                         chr = seqnames(probeInfo)[hits@to],
                         CpG = start(probeInfo)[hits@to],
                         CGIposition = probeInfo@elementMetadata$CGIposition[hits@to],
                         stringsAsFactors = F)



head(ROItoProbe)

######

geneToProbe <- merge(geneToROI,ROItoProbe,by="ROI",all.x=T,sort=F)

inshores <- which(geneToProbe$CGIposition %in% c("N_Shore","S_Shore"))
geneToProbe$ROI[inshores] <- gsub("CGI","SHO",geneToProbe$ROI[inshores])

geneToProbe$distToTSS <- geneToProbe$CpG - geneToProbe$transcription_start_site

apply(is.na(geneToProbe),2,sum)

SHO <- ROI[substr(names(ROI),1,3) == "CGI"]
names(SHO) <- gsub("CGI","SHO",names(SHO))

ROI <- c(ROI,SHO)
ROI <- sort(ROI)

## ---- auxiliary functions ----

plotROI <- function(roi,viewportadd=5e3,geneNames="",...) {
  x <- ROI[roi]
  
  viewport <- x
  start(viewport) <- start(x) - viewportadd
  end(viewport) <- end(x) + viewportadd
  
  ensts <- ENST[findOverlaps(ENST,viewport,ignore.strand=T)@from]
  genes <- ensts$hgnc_symbol
  cgis <- CGI[findOverlaps(CGI,viewport,ignore.strand=T)@from]
  probes <- probeInfo[findOverlaps(probeInfo,viewport,ignore.strand=T)@from]
  
  plot(NA,xlim=c(start(viewport),end(viewport)),ylim=c(6.5,0),
       xlab=sprintf("%s coordinate",seqnames(x)),yaxt="n",ylab="",...)
  
  rect(start(x)-2,1,end(x)+2,1.5,col="lightblue",border=NA)
  text(mean(c(start(x),end(x))),.5,roi,cex=1)
  if(length(cgis) > 0) rect(start(cgis)-2,1.9,end(cgis)+2,2.7,col="pink",border=NA)
  
  segments(start(probes),2.1,start(probes),2.5)
  
  if(length(ensts) > 0) {
    ymax <- min(c(6,3+length(ensts)*.2)) # space the ENSTs between tracks 3 and 6
    colors <- ifelse(genes %in% geneNames,"black","grey")
    widths <- ifelse(genes %in% geneNames,2,1)
    yensts <- seq(3,ymax,l=length(ensts))
    segments(start(ensts),yensts,end(ensts),yensts,col=colors,lwd=widths)
    tss <- ifelse(ensts@strand=="+",start(ensts),end(ensts))
    tts <- ifelse(ensts@strand=="+",end(ensts),start(ensts))
    
    points(tss,yensts,pch=21,bg="green")
    points(tts,yensts,pch=21,bg="red")
  }
  
  
}



## save data and functions ----

save(PROMOTERWIDTH,ENST,ROI,CGI,pbTSS,geneToProbe,plotROI,file="data/regions.rda")
