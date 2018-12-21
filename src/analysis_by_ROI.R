# Classify samples and ROIs
# By Sergio Alonso

library(gplots)
library(corrplot)
library(dplyr)
library(igraph)
library(RColorBrewer)

cancer <- "CRC" # testing purposes. This must be commented/removed in the final version

methylationdir <- "/imppc/labs/mplab/share/ICGC/"
wd <- "/imppc/labs/mplab/share/metallopeptidases/"

setwd(wd)
somatic <- readRDS(sprintf("%s/%s/methICGC/somaticROI.rds",methylationdir,cancer))
load("data/pantherdb.rda",v=T)
load("data/regions.rda",v=T)
load("data/probeInfo.rda",v=T)
distances <- readRDS("data/distances.rds")
source("src/checkBimodality.R")

analyzeSomatic <- function(somatic,genes=paste0("ADAMTS",c(1:10,12:19)),
                           geneColors = NULL,
                           ROITypes=c("CGI","SHO","TSS"),
                           cisThreshold = 1e6,
                           selectionThreshold=0.05,
                           generatePlots=T,
                           tumorOrder = c("local,globalCGI")) {
  
 
  # rois is a dataframe with the information for the analyzed ROIs
  
  rois <- subset(geneToProbe,hgnc_symbol %in% genes)
  rois <- rois[,c("hgnc_symbol","ROI")]
  rois$type <- substr(rois$ROI,1,3)
  rois <- subset(rois,type %in% ROITypes) %>% unique
  
  if(is.null(geneColors)) {
    rois$color <- factor(rois$type,c("CGI","SHO","TSS"),c("red","orange","blue"))
    } else {
    names(geneColors) <- genes
    rois$color <- geneColors[rois$hgnc_symbol]
  }
  
  rois <- rois[rois$ROI %in% rownames(somatic),] # only ROIs with info in the somatic data table

  genenames <- tapply(rois$hgnc_symbol,rois$ROI,function(x) x[x!=""] %>% paste(.,collapse=";"))
  
  rois$genenames <- genenames[rois$ROI]
  
  rm(genenames)
  
  namelength <- nchar(rois$ROI) + nchar(rois$genenames)

  spaces <- sapply(namelength,function(x) paste(rep(" ",max(namelength)-x+1),collapse=""))

  rois$longnames <- sprintf("%s%s%s",rois$ROI,spaces,rois$genenames)
  
  rois$hgnc_symbol <- NULL
  
  newcolors <- tapply(rois$color,rois$ROI,function(x) {
    color <- unique(x)
    color <- colorRamp(color)(0.5)
    rgb(color[1],color[2],color[3],maxColorValue = 255)
  }) 
  
  rois$color <- newcolors[rois$ROI]
  
  rois <- unique(rois)
  rownames(rois) <- rois$ROI
  
  # prepare the datatable with the somatic data of the gene-associated ROIs

  tumors <- substring(colnames(somatic),14,15) == "01"
  somatic0 <- somatic[rois$ROI,tumors]
  
  cgimeth <- colMeans(somatic[substr(rownames(somatic),1,3)=="CGI",tumors],na.rm=T) # meytylation in all CGIs
  roimeth <- colMeans(somatic0,na.rm=T) # methylation in ROIs

  corGlobal <- cor(cgimeth,t(somatic0),use="complete")
  
  switch(match.arg(tumorOrder,c("local","globalCGI")),
         
         local = ordertumors <- order(colMeans(somatic0,na.rm=T)),
         globalCGI = ordertumors <- order(cgimeth) 
         
         )
         
  orderROIs <- order(rowMeans(somatic0,na.rm=T))

  somaticChanges <- apply(abs(somatic0) > .25,1,sum,na.rm=T) 
  selected <- somaticChanges >= ncol(somatic) * selectionThreshold # altered in 5% of the samples 
  
  # x: somatic alterations selected
  # cor.x: corraltion of x
  # cor.x2: correlation of x, ordered according to hclust complete
  # and with in-cis correlations eliminated
  
  x <- somatic0[selected,]
  cor.x <- cor(t(x),use = "complete")
  
  h1 <- hclust(as.dist(1 - cor.x),method = "complete")
  d1 <- as.dendrogram(h1)
  clusters <- cutree(h1,h = 1-.25)
  
  order1 <- h1$order
  
  cor.x2 <- cor.x[order1,order1]
  
  cis <- cor.x2 # to copy size and row/colnames
  cis[,] <- cisThreshold+10 # set all values to a value larger than cisThreshold
  
  for(i in rownames(cis)) {
    chr <- as.character(seqnames(ROI[i]))
    d <- distances[[chr]][i,]
    inChr <- intersect(colnames(cis),names(d))
    cis[i,inChr] <- d[inChr]
  }
  
 
  cor.x2[cis < cisThreshold] <- NA
  
  rownames(cor.x2) <- rois[rownames(cor.x2),"longnames"]
  
  
  if(generatePlots) {
    par(family="mono",font=1)
    
    
    # plot the analyzed ROIs
    
    
    if(nrow(rois) <= 4) rowcol <- c(nrow(rois),1)
    if(nrow(rois) > 4) rowcol <- c(4,2)
    
    par(mfrow=rowcol,mar=c(4.5,4,1,.5),mgp=c(2.3,.8,0),xpd=F)
    
    for(i in 1:nrow(rois)) plotROI(roi = rois$ROI[i],
                                   geneName = unlist(strsplit(rois$genenames[i],split = ";")),
                                   main=rois$genenames[i])
    
    # Heatmap with all analyzed ROIs
    
    par(mar=c(5,5,2,1),mfrow=c(1,1))
    
    somaticPalette <- colorRampPalette(c("darkblue","white","darkred"))(20)
    correlationPalette <- colorRampPalette(c("orange4","white","darkgreen"))(20)
    somaticRamp <- colorRamp(c("darkblue","white","darkred"))
    correlationRamp <- colorRamp(c("orange4","white","darkgreen"))
  
    tumorColor <- somaticRamp(cgimeth/max(abs(cgimeth))/2 + .5)/255 
    tumorColor <- apply(tumorColor,1,function(x) rgb(x[1],x[2],x[3]))
    tumorColor <- tumorColor
    
    roiColor <- correlationRamp(corGlobal/2 + .5) / 255
    roiColor <- apply(roiColor,1,function(x) rgb(x[1],x[2],x[3]))
    
    heatmap.2(somatic0[orderROIs,ordertumors],trace="none",dendrogram = "none",
              Colv = NULL,
              Rowv = NA,
              margins = c(5,30),
              breaks=seq(-.9,.9,l=21),
              col = somaticPalette,
              labRow = rois$longnames[orderROIs],
              colRow = NA,
              labCol=NA,
              keysize=1,
              key.title = "Somatic alterations",
              density.info="hist",
              cexRow = 1,
              cexCol = 1,
              ColSideColors = tumorColor[ordertumors],
              RowSideColors = roiColor[orderROIs])
  
  
    heatmap.2(somatic0[orderROIs,ordertumors],trace="none",dendrogram = "none",
              Colv = NULL,
              Rowv = NA,
              margins = c(5,30),
              breaks=c(-1,-.25,.25,1),
              col = colorRampPalette(c("darkblue","white","darkred"))(3),
              labRow = rois$longnames[orderROIs],
              colRow = NA,
              labCol=NA,
              keysize=1,
              key.title = "Somatic alterations",
              density.info="hist",
              cexRow = 1,
              cexCol = 1,
              ColSideColors = tumorColor[ordertumors],
              RowSideColors = roiColor[orderROIs])
    
  
  # Heatmap with selected ROIs
  
    nclusters <- nlevels(factor(clusters))
    clusterColors <- brewer.pal(nclusters,"Set3")
    
      heatmap.2(x[,ordertumors],trace="none",dendrogram = "row",
              Colv = NA,
              Rowv = rev(d1),
              margins = c(5,30),
              breaks = seq(-.9,.9,l=21),
              col = somaticPalette,
              labRow = rois$longnames[selected],
              labCol = NA,
              keysize = 1.2,
              key.title = "Somatic alterations",
              cexRow = 1,
              cexCol = 1,
              ColSideColors = tumorColor[ordertumors],
              RowSideColors = clusterColors[clusters])
      
      heatmap.2(x[,ordertumors],trace="none",dendrogram = "row",
                Colv = NA,
                Rowv = rev(d1),
                margins = c(5,30),
                breaks = c(-1,-.25,.25,1),
                col = colorRampPalette(c("darkblue","white","darkred"))(3),
                labRow = rois$longnames[selected],
                labCol = NA,
                keysize = 1.2,
                key.title = "Somatic alterations",
                cexRow = 1,
                cexCol = 1,
                ColSideColors = tumorColor[ordertumors],
                RowSideColors = clusterColors[clusters])
  
  
  # Heatmap of correlations
  
    heatmap.2(cor.x,trace="none",dendrogram = "both",
              Colv = d1,
              Rowv = rev(d1),
              margins = c(5,30),
              breaks = seq(-1,1,l=21),
              col = correlationPalette,
              labRow = rois$longnames[selected],
              labCol = NA,
              key.title = "Correlation",
              keysize=1.2,
              density.info="h",
              cexRow = 1,
              cexCol = 1,
              ColSideColors = clusterColors[clusters],
              RowSideColors = clusterColors[clusters])
 
  
  
  # Corrplot of correlations (deleting correlations in cis)
  # find in cis correlations
  
    corrplot(cor.x2,method="square",type="full",diag=T,col=correlationPalette,
             tl.col="black",
             na.label.col="black",na.label = ".",
             mar=c(2,2,2,2))
    
    corrRect.hclust(cor.x,k=nclusters,method = "complete")
  
  
  # igraph plot of the correlations
  
    ig0 <- graph_from_adjacency_matrix(cor.x2,mode="upper",weighted = T,diag=F)
    E(ig0)$correlation <- E(ig0)$weight
    ig0 <- delete_edges(ig0,which(is.na(E(ig0)$correlation)))
    ig0 <- delete_edges(ig0,which(abs(E(ig0)$correlation) < 0.25))
    
    E(ig0)$weight <- 1+E(ig0)$weight
    E(ig0)$color <- apply(correlationRamp((E(ig0)$correlation+1)/2),1,function(x) rgb(x[1],x[2],x[3],maxColorValue = 255))
    
    E(ig0)$width <- abs(E(ig0)$correlation) * 10
    E(ig0)$curved <- .25
    V(ig0)$color <- as.character(rois[V(ig0)$name,"color"])
    V(ig0)$label <- rois[V(ig0)$name,"genenames"]
    V(ig0)$label.color <- "black"
    V(ig0)$label.cex <- .5
    V(ig0)$cex <- .5
    communities0 <- cluster_label_prop(ig0)
    plot.igraph(ig0,layout=layout_with_fr,mark.groups = communities0,mark.border=NA,vertex.label=NA)

    ig1 <- delete_edges(ig0,which(abs(E(ig0)$correlation) < 0.5))
    ig1 <- delete_vertices(ig1,which(degree(ig1)==0))
    communities1 <- cluster_label_prop(ig1)
    plot.igraph(ig1,layout=layout_with_fr,mark.groups = communities1,mark.border=NA,vertex.label=NA)
    
    
  }
  
  clusters2 <- list()
  
  for(i in 1:max(clusters)) {
    inCluster <- which(colnames(cor.x2) %in% names(clusters)[clusters==i])
    clusters2[[i]] <- cor.x2[inCluster,inCluster] %>% .[upper.tri(.)] %>% .[!is.na(.)]
  } 
  
  names(clusters2) <- paste0("Cluster",1:max(clusters))
  
  granularity <- length(unlist(clusters2)) / (length(clusters) * (length(clusters) - 1) / 2)
  
  invisible(list(genes=genes,
                 ROIs=rois$ROI,
                 somatic=x,
                 correlations=cor.x2,
                 clusters=clusters,
                 clustersCorr=clusters2,
                 Gr = granularity,
                 mGr = mean(unlist(clusters2))))
  
}


genes <- sapply(c("ADAMTS1","MMP1","ADAM8"),function(x) familyToGene(geneToFamily(x)))
colors <- rep(c("pink","lightblue","lightgreen"),sapply(genes,length))
genes <- as.character(unlist(genes))

pdf("sandbox/metallopeptidases.pdf",width = 10,height = 10)

foo <- analyzeSomatic(somatic,
                      genes = genes,
                      #geneColors = colors,
                      ROITypes = c("CGI","TSS"),
                      selectionThreshold = 0.05,
                      generatePlots = T,
                      tumorOrder = "global")

dev.off()


familyToGene(geneToFamily("MMP1"))




# select genes and analysis to perform

analyses <- list()

analyses$ADAMTS <- familyToGene(geneToFamily("ADAMTS19"))
analyses$ADAM <- familyToGene(geneToFamily("ADAM2"))
analyses$MMP <- familyToGene(geneToFamily("MMP8"))
analyses$HOX <- familyToGene(geneToFamily("HOXA1"))

x <- unique(subset(panther2,Freq > 10 & Freq < 50)$hmmpanther)
y <- lapply(x,familyToGene)
names(y) <- x

analyses <- c(analyses,y)

lapply(analyses,function(analysis) {
  foo <- NA
  try(foo <- analyzeSomatic(somatic=somatic,
                            genes=analysis,
                            ROITypes = c("CGI","TSS"),
                            generatePlots = F,
                            cisThreshold = 10e6,
                            selectionThreshold = 0.05,
                            tumorOrder = "global"))
  foo
}) -> results

# remove results that failed

results <- results[sapply(results,function(x) !is.null(names(x)))]

library(Biobase)


final <- data.frame(genes = sapply(subListExtract(results,"genes"),length),
                    ROIs = sapply(subListExtract(results,"ROIs"),length),
                    somROIs = sapply(subListExtract(results,"clusters"),length),
                    cis = sapply(subListExtract(results,"correlations"),function(x) {
                      x <- x[upper.tri(x)]
                      sum(is.na(x)) / length(x)
                    }),
                    mCorr = sapply(subListExtract(results,"correlations"),mean,na.rm=T),
                    Gr = subListExtract(results,"Gr",simplify = T),
                    mGr = subListExtract(results,"mGr",simplify = T))

niceplot <- function(x,y,selected=1:4,colors=c("red","blue","green","orange"),...) {
  plot(x,y,pch=19,col="lightblue",...)
  grid()
  points(x[selected],y[selected],pch=21,bg=colors,cex=1.5)
  text(x[selected],y[selected],rownames(final)[selected],pos=4)
  sfam <- identify(x,y,labels = panther2[rownames(final),"name"],cex=.8)
  points(x[sfam],y[sfam],pch=1,cex=1.2)
  print(panther2[rownames(final)[sfam],])
  }

niceplot(final$genes,final$ROIs,xlab="Genes in family",ylab="ROIs",xlim=c(10,60),las=1)
niceplot(final$ROIs,final$somROIs/final$ROIs*100,xlab="ROIs",ylab="ROIs with alterations (%)",las=1)
niceplot(final$somROIs,final$cis*100,xlab="ROIs with somatic alterations",ylab="cis interactions (%)",las=1)
niceplot(final$somROIs/final$ROIs*100,final$cis*100,xlab="ROIs with alterations (%)",ylab="cis interactions (%)",las=1)


niceplot(final$somROI,final$Gr,xlab="ROIs with somatic alterations",ylab="Granularity")
niceplot(mGr ~ somROI,xlab="ROIs with somatic alterations",ylab="In-cluster correlation",ylim=c(0,1))


