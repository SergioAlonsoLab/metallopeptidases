library(shiny)
library(dplyr)
library(igraph)
library(RColorBrewer)
library(karyoploteR)
library(gplots)
library(corrplot)

setwd("/imppc/labs/mplab/share/metallopeptidases/")
load("data/pantherdb.rda",v=T)
load("data/regions.rda",v=T)
load("data/probeInfo.rda",v=T)

cancers <- c("CRC","BRCA","LUAD")

geneToProbe$ROIType <- factor(substr(geneToProbe$ROI,1,3),c("TSS","CGI","SHO"))
ROIcolors <- c("lightblue","red","orange")

x <- sapply(unique(geneToProbe$ROI), function(x) strsplit(x,split="(-|:)")[[1]][3:4] %>% as.numeric) %>% t
geneToProbe$ROIstart <- x[geneToProbe$ROI,1]
geneToProbe$ROIend <- x[geneToProbe$ROI,2]
rm(x)

famCol <- rep(RColorBrewer::brewer.pal(9,"Pastel1"),10)
rfl <- removeFirstLast <- function(x) x[-c(1,length(x))]


# SHINY LAYOUT ----

ui <- fluidPage(
  tabsetPanel(id = "tabs",
              
              tabPanel("Cancer Type", 
                       selectInput("cancerType","Cancer Type",choices = cancers,selected = "CRC")),
              
              tabPanel("Genes and Families", 
                       fluidRow(
                         
                         column(width = 3,
                                textAreaInput("genes","Genes",height = "500px",resize="none"),
                                actionButton("toFamilies","To Families"),
                                actionButton("clearGenes","Clear")
                         ),
                         
                         column(width = 3,
                                textAreaInput("families","Gene Family",height = "500px",resize="none"),
                                actionButton("toGenes","To Genes"),
                                actionButton("removeSF","Remove SubFamilies"),
                                actionButton("clearFamilies","Clear")
                         ),
                         
                         column(width = 6,
                                textAreaInput("familyNames","Family Name",height = "500px",resize="none"))
                       )
              ),
              
              
              tabPanel("Genes on Karyotype",
                       plotOutput("karyotype",width="1200px",height="900px")),
              
              tabPanel("ROIs", 
                       sidebarLayout(
                         sidebarPanel(width = 3,
                                      checkboxGroupInput("ROITypes",label = "ROI Type",
                                                         choiceNames = list("TSS","CGI","Shores"),
                                                         choiceValues = list("TSS","CGI","SHO"),
                                                         selected = list("TSS","CGI")),
                                      
                                      actionButton("updateROIs","Update ROIS")),
                         
                         mainPanel(width = 9,
                                   plotOutput("grois",height="800px",width="1000px"))
                         
                       )),
              
              tabPanel("ROIs Location", 
                       sidebarLayout(
                         
                         sidebarPanel(width = 3,
                                      selectInput("selectedGene","GENE",choices = "NONE",multiple = F)),
                         
                         mainPanel(width = 9,
                                   plotOutput("roisPlot",width = "1000px",height = "400px"))
                         
                         
                       )
              ),
              
              
              tabPanel("Somatic alterations",
                       sidebarLayout(
                         sidebarPanel(width = 2,
                                      sliderInput("somaticBins","Number of Bins",3,51,15,step = 2),
                                      sliderInput("selThreshold","Selection Threshold",0,.8,0,step=.05)),
                         
                         mainPanel(width = 10,
                                   plotOutput("gsom",width="1400px",height="1000px"))
                         
                       )),
              
              tabPanel("Correlation graph A",
                       sidebarLayout(
                         sidebarPanel(width = 2,
                                      numericInput("cisThreshold","cis Threshold",value = 1e6,min = 0, max = 1e9, step = 1000)),
                         mainPanel(width = 10,
                                   plotOutput("corplot1",width="1400px",height="1000px"))
                       )),
              
              
              
              
              
              
              
              tabPanel("Correlation graph B",
                       sidebarLayout(
                         sidebarPanel(width = 2,
                                      sliderInput("clusterThreshold","Min corr cluster",value=0.25, min = 0,max = 1, step = 0.05),
                                      textOutput("granularity")),
                         mainPanel(width = 10,
                                   plotOutput("corplot2",width="1400px",height="1000px"))
                       ))
              
              
  ))



# SHINY SERVER -------


server <- function(input,output,session) {
  
  # color palettes
  
  somaticPalette <- colorRampPalette(c("darkblue","white","darkred"))
  correlationPalette <- colorRampPalette(c("orange4","white","darkgreen"))
  somaticRamp <- colorRamp(c("darkblue","white","darkred"))
  correlationRamp <- colorRamp(c("orange4","white","darkgreen"))
  
  
  
  # select the cancer type
  
  methylationdir <- "/imppc/labs/mplab/share/ICGC/"
  cancerData <- reactive({
    cancer <- input$cancerType
    somatic <- readRDS(sprintf("%s/%s/methICGC/somaticROI.rds",methylationdir,cancer))
    tumors <- substring(colnames(somatic),14,15) == "01"
    cgimeth <- colMeans(somatic[substr(rownames(somatic),1,3)=="CGI",tumors],na.rm=T)
    cat("updated cancer data\n")
    return(list(somatic=somatic,tumors=tumors,cgimeth=cgimeth))
  })
  
  
  # rois is the main data frame where the genes, mirnas, probes, etc are stored
  
  rois <- reactive({
    genes <- readInput(input$genes)
    ROITypes <- input$ROITypes
    
    rois <- subset(geneToProbe,hgnc_symbol %in% genes & ROIType %in% ROITypes) # this is the main information table
    nameLength <- nchar(rois$ROI) + nchar(rois$hgnc_symbol)
    rois$longName <- sprintf("%s%s%s",rois$ROI,strrep(".",max(nameLength)-nameLength+2),rois$hgnc_symbol)
    
    cat("updated rois\n")
    return(rois)
    
  })
  
  
  somatic <- reactive({
    cancerData <- cancerData()
    rois <- rois()
    somatic <- cancerData$somatic[intersect(rownames(cancerData$somatic),rois$ROI),cancerData$tumors]
    
    selected <- rowSums(abs(somatic) > input$selThreshold) >= sum(tumors)*.05
    cat("updated somatic data\n")
    return(somatic[selected,])
    
  })
  
  
  correlations <- reactive({
    
    somatic <- somatic()
    rois <- rois()
    cisThreshold <- input$cisThreshold
    
    longNames <- rois$longName
    names(longNames) <- rois$ROI
    
    cor.x2 <- cor.x <- cor(t(somatic),use="complete")
    
    h1 <- hclust(as.dist(1 - cor.x),method = "complete")
    d1 <- as.dendrogram(h1)
    clusters <- cutree(h1,h = 1-.25)
    
    selectedROIs <- ROI[rownames(somatic)]
    d <- sapply(1:length(selectedROIs),function(x) distance(selectedROIs[x],selectedROIs))
    d[is.na(d)] <- cisThreshold + 1
    
    cor.x2[d < cisThreshold] <- NA
    
    cat("updated correlations\n")
    
    return(list(cor.x=cor.x,
                cor.x2=cor.x2,
                h1=h1,
                d1=d1,
                clusters=clusters,
                longNames=longNames))
    
    
  })
  
  
  granularity <- reactive({
    
    with(correlations(),{
      clusters <- cutree(h1,h = 1-input$clusterThreshold)
      
      nclusters <- nlevels(factor(clusters))
      valid <- sum(!is.na(cor.x2))
      inCluster <- sapply(1:nclusters,function(i) sum(!is.na(cor.x2[clusters==i,clusters==i]))) %>% sum
      return(inCluster/valid)
    })
    
    
  })
  
  
  
  observeEvent(input$toGenes,{
    families <- readInput(input$families)
    genes <- familyToGene(families)
    updateGenesAndFamilies(genes,families,session)
  })
  
  observeEvent(input$removeSF, {
    families <- readInput(input$families)
    families <- families[!grepl(":SF",families)]
    genes <- readInput(input$genes)
    updateGenesAndFamilies(genes,families,session)
  })
  
  observeEvent(input$clearFamilies,{
    updateTextAreaInput(session,"families",value="")
    updateTextAreaInput(session,"familyNames",value="")
  })
  
  observeEvent(input$toFamilies,{
    genes <- readInput(input$genes)
    families <- geneToFamily(genes)
    updateGenesAndFamilies(genes,families,session)
  })
  
  observeEvent(input$clearGenes,{
    updateTextAreaInput(session,"genes",value="")
  })
  
  
  observeEvent(input$updateROIs, {
    families <- readInput(input$families)
    nmax <- sapply(list(families,genes,rois()$ROI), function(x) length(unique(x))) %>% max
    output$grois <- renderPlot(
      height = 800+(nmax-20)*10,
      expr = roisGraph(rois(),families) 
    )
  })
  
  
  # plot the graphs
  
  output$karyotype <- renderPlot(mapgenes(readInput(input$genes)))
  output$roisPlot <- renderPlot(mapROIs(input$selectedGene,rois()))
  output$gsom <- renderPlot(plotSomatic(somatic(),cancerData()$cgimeth,rois(),input$somaticBins))
  output$corplot1 <- renderPlot(plotCorrelations(correlations()))
  output$corplot2 <- renderPlot(plotCorrelations2(correlations(),input$clusterThreshold))
  output$granularity <- renderText(sprintf("In clusters: %1.2f%%",granularity()*100))
  
}


# Functions to be run on the Shiny Server ----

readInput <- function(inputId) {
  x <- strsplit(inputId,split="(\n|,)")[[1]] %>% gsub(" +","",.)
  x[x!=""]
}

updateGenesAndFamilies <- function(genes,families,session) {
  familyNames <- paste(families,panther2[families,"name"],sep=": ")
  updateTextAreaInput(session,"genes",value=paste(genes,collapse="\n"))
  updateTextAreaInput(session,"families",value=paste(families,collapse="\n"))
  updateTextAreaInput(session,"familyNames",value=paste(familyNames,collapse="\n\n"))
  updateSelectInput(session = session,inputId = "selectedGene",choices=genes)
}

updateROIs <- function(genes,ROITypes) {
  rois <- subset(geneToProbe,hgnc_symbol %in% genes & ROIType %in% ROITypes) # this is the main information table
  nameLength <- nchar(rois$ROI) + nchar(rois$hgnc_symbol)
  rois$longName <- sprintf("%s%s%s",rois$ROI,strrep(".",max(nameLength)-nameLength+2),rois$hgnc_symbol)
  return(rois)
}

# Graph showing the families, genes and ROIs ----

roisGraph <- function(rois,families) {
  
  
  print(colnames(rois))
  rois <- rois[,c("hgnc_symbol","ensembl_transcript_id","ROI","chr","longName")]
  rois <- unique(rois)
  
  # add columns with family information
  # one column per family, with T or F
  
  m2 <- sapply(families,function(fam) {
    sapply(unique(rois$hgnc_symbol),function(x) fam %in% geneToFamily(x))
  })
  
  m2 <- as.matrix(m2)
  if(ncol(m2) > 1) m2 <- m2[,order(colSums(m2))] # order them from less to more members
  
  # Dirty trick to order genes according to its family
  # there is probably a more elegant way, but this works
  o1 <- do.call("order",as.data.frame(m2[,ncol(m2):1]))
  o2 <- 1:length(o1)
  names(o2) <- rownames(m2)[rev(o1)]
  
  rois$order <- o2[rois$hgnc_symbol]
  rois <- rois[order(rois$order,rois$hgnc_symbol),] # order ROIs according to the gene family
  
  m <- c(families=list(rownames(m2)),apply(rois,2,unique))
  m$type <- factor(substr(m$ROI,1,3),c("TSS","CGI","SHO"))
  n <- sapply(m,length)
  
  plot.new()
  
  x0 <- min(0,(10-ncol(m2))/10)
  
  plot.window(xlim=c(x0,5),ylim=c(10,-1))
  
  rfl <- removeFirstLast <- function(x) x[-c(1,length(x))]
  
  y <- lapply(n,function(x) rfl(seq(0,10,l=x+2)))
  
  w <- c(.5,.8,.2,.4)/2
  for(i in 2:4) {
    h <- min(1,10/n[i])/3
    col <- NA
    if(i == 3) col <- "lightgreen"
    if(i == 4) col <- c("lightblue","red","orange")[m$type] 
    rect(i-w[i],y[[i]]-h,i+w[i],y[[i]]+h,col=col)
  }
  
  text(2:4,-.2,c("Genes","mRNAs","ROIs"))
  text(2,y[[2]],m$hgnc_symbol)
  
  x <- 2-w[2]-(1:(ncol(m2)))/10
  h <- min(1,10/n[2])/3
  
  for(i in 1:ncol(m2)) {
    sel <- sapply(m$hgnc_symbol,function(x) colnames(m2)[i] %in% geneToFamily(x))
    rect(x[i],y[[2]]-h,x[i]-.1,y[[2]]+h,col=ifelse(sel==T,famCol[i],"white"))
  }
  
  text(x-.05,0,colnames(m2),srt=90,adj=0,cex=1)
  
  for(i in 1:n[3]) {
    sel <- m$hgnc_symbol %in% rois$hgnc_symbol[rois$ensembl_transcript_id==m$ensembl_transcript_id[i]]
    segments(2+w[2],y[[2]][sel],3-w[3],y[[3]][i])
  }
  
  for(i in 1:n[4]) {
    sel <- m$ensembl_transcript_id %in% rois$ensembl_transcript_id[rois$ROI==m$ROI[i]]
    segments(3+w[3],y[[3]][sel],4-w[4],y[[4]][i])
  }
  
  
}

# map genes on Karyotype

mapgenes <- function(genes) {
  enst <- subset(ENST,hgnc_symbol %in% genes)
  tss <- ifelse(strand(enst)=="+",start(enst),end(enst))
  minTSS <- tapply(tss,enst$hgnc_symbol,min)
  maxTSS <- tapply(tss,enst$hgnc_symbol,max)
  chr <- tapply(as.character(seqnames(enst)),enst$hgnc_symbol,unique)
  
  kp <- plotKaryotype()
  if(length(genes) > 0) {
    kpPlotMarkers(kp,chr=chr,x=(minTSS+maxTSS)/2,labels=names(minTSS),text.orientation = "hor",y=.5,cex=.8)
  }
  
}


# Individual plots of the loctaion of the ROIs

mapROIs <- function(gene,rois) {
  x <- subset(rois,hgnc_symbol == gene)
  ensts <- ENST[unique(x$ensembl_transcript_id)]
  region <- range(ensts) + 5000
  geneStartEnd <- range(c(start(ensts),end(ensts)))
  geneLength <- diff(geneStartEnd)
  strand <- ifelse(unique(strand(ensts))=="+",1,-1)
  geneStart <- ifelse(strand == 1,geneStartEnd[1],geneStartEnd[2])
  
  
  plot(NA,xlim=c(start(region),end(region)),ylim=c(0.5,3.5),
       yaxt="n",
       xlab=sprintf("%s:%i-%i",unique(x$chr),start(region),end(region)),
       ylab="")
  rect(geneStart,1,geneStart + strand * geneLength*.91,1.2,border=NA,col="lightblue")
  polygon(geneStart + strand * geneLength * c(.9,1,.9),c(.95,1.1,1.25),border=NA,col="lightblue")
  
  text(mean(geneStartEnd),1,sprintf(ifelse(geneLength < 1e6,sprintf("%1.1f Kb",geneLength/1000),sprintf("%1.1f Gb",geneLength/1e6))),pos=1)
  text(mean(geneStartEnd),1.1,gene,font=3)
  y.ensts <- removeFirstLast(seq(1.2,2.5,l=length(ensts)+2))
  
  ensts.start <- ifelse(strand(ensts)=="+",start(ensts),end(ensts))
  ensts.stop <- ifelse(strand(ensts)=="+",end(ensts),start(ensts))
  
  segments(ensts.start-PROMOTERWIDTH/2,y.ensts,ensts.start+PROMOTERWIDTH/2,y.ensts,lwd=6,col="grey")
  segments(ensts.start,y.ensts,ensts.stop,y.ensts,lty=ifelse(ensts$transcript_tsl=="tsl1",1,2))
  points(ensts.start,y.ensts,pch=21,bg="green")
  points(ensts.stop,y.ensts,pch=21,bg="red")
  
  cgis <- CGI[findOverlaps(CGI,region)@from]
  if(length(cgis) > 0) rect(start(cgis),2.5,end(cgis),2.7,col="pink",border="pink")
  
  probes <- probeInfo[findOverlaps(probeInfo,region,ignore.strand=T)@from]
  segments(start(probes),2.8,end(probes),2.9)
  
  rect(start(region)-6000,3-.05,end(region)+6000,3.3+.05,col="lightgrey",border=NA)
  y <- unique(x[,c("ROI","ROIType")])
  rois0 <- ROI[as.character(unique(x$ROI))]
  rect(start(rois0),3,end(rois0),3.3,col=ROIcolors[y$ROIType])
  
  
}


#### SOMATIC CHANGES #####

plotSomatic <- function(somatic,cgimeth,rois,nbins) {
  
  roimeth <- colMeans(somatic,na.rm=T) # methylation in ROIs
  corGlobal <- cor(cgimeth,t(somatic),use="complete")
  orderROIs <-order(rowMeans(somatic,na.rm=T))
  ordertumors <- order(cgimeth) 
  
  par(family="mono",font=1)
  
  tumorColor <- somaticRamp(cgimeth/max(abs(cgimeth))/2 + .5)/255 
  tumorColor <- apply(tumorColor,1,function(x) rgb(x[1],x[2],x[3]))
  tumorColor <- tumorColor
  
  roiColor <- correlationRamp(corGlobal/2 + .5) / 255
  roiColor <- apply(roiColor,1,function(x) rgb(x[1],x[2],x[3]))
  
  longNames <- rois$longName
  names(longNames) <- rois$ROI
  labelcex <- min(.1 + 100/length(orderROIs),3)
  
  heatmap.2(somatic[orderROIs,ordertumors],trace="none",dendrogram = "none",
            Colv = NULL,
            Rowv = NA,
            margins = c(5,60),
            breaks=seq(-.9,.9,l=nbins+1),
            col = somaticPalette(nbins),
            labRow = longNames[rownames(somatic)][orderROIs],
            colRow = NA,
            labCol=NA,
            keysize=1,
            key.title = "Somatic alterations",
            density.info="hist",
            cexRow = max(2.5,.2 + 2.5/log10(nrow(somatic))),
            cexCol = 1,
            ColSideColors = tumorColor[ordertumors],
            RowSideColors = roiColor[orderROIs])
  
}

# Correlations graphs ---------

plotCorrelations <- function(correlations) {
  par(family="mono",font=1)
  
  with(correlations,{
    
    nclusters <- nlevels(factor(clusters))
    clusterColors <- rainbow(n = nclusters,s = .3, v = 1)
    
    heatmap.2(cor.x2,trace="none",dendrogram = "both",
              Colv = d1,
              Rowv = rev(d1),
              margins = c(5,60),
              #breaks = seq(-1,1,l=21),
              col = correlationPalette,
              labRow = longNames[rownames(cor.x)],
              labCol = NA,
              key.title = "Correlation",
              keysize=1.2,
              density.info="h",
              cexRow = max(2.5,.2 + 2.5/log10(nrow(cor.x2))),
              #cexCol = 1,
              ColSideColors = clusterColors[clusters],
              RowSideColors = clusterColors[clusters],
              na.color="grey50")
  })
  
}


plotCorrelations2 <- function(correlations,clusterThreshold) {
  par(family="mono",font=1)
  
  with(correlations,{
    clusters <- cutree(h1,h = 1-clusterThreshold)
    
    nclusters <- nlevels(factor(clusters))
    clusterColors <- rainbow(n = nclusters,s = .3, v = 1)
    
    o1 <- h1$order
    rownames(cor.x2) <- longNames[rownames(cor.x2)]
    corrplot(cor.x2[o1,o1],method="square",type="full",diag=T,col=correlationPalette,
             tl.col="black",
             na.label.col="black",na.label = "X",
             mar=c(2,2,2,2))
    
    corrRect.hclust(cor.x,k=nclusters,method = "complete")
    
    
    
  })
}

# granularity value



shinyApp(ui,server)
