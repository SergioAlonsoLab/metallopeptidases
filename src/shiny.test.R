library(shiny)
library(dplyr)
library(igraph)
library(RColorBrewer)
library(karyoploteR)

setwd("/imppc/labs/mplab/share/metallopeptidases/")
load("data/pantherdb.rda",verbose = T)
load("data/regions.rda",v=T)

geneToProbe$ROIType <- substr(geneToProbe$ROI,1,3)
famCol <- rep(RColorBrewer::brewer.pal(9,"Pastel1"),10)


# SHINY LAYOUT ----

ui <- fluidPage(
  tabsetPanel(
    
    tabPanel("Cancer Type", 
             selectInput("cancerType","Cancer Type",choices = cancers,selected = "CRC")),
    
    tabPanel("Genes and Families", {
      fluidRow(
        
        column(width = 3,
               textAreaInput("genes","Genes",height = "500px",resize="none"),
               column(width = 6, actionButton("toFamilies","To Families")),
               column(width = 6, actionButton("clearGenes","Clear"))
        ),
        
        column(width = 3,
               textAreaInput("families","Gene Family",height = "500px",resize="none"),
               column(width = 6,actionButton("toGenes","To Genes")),
               column(width = 6,actionButton("clearFamilies","Clear"))
        ),
        
        column(width = 6,
               textAreaInput("familyNames","Family Name",height = "500px",resize="none"))
      )
    }),
    
    tabPanel("ROIS", {
      sidebarLayout(
        sidebarPanel(width = 3,
               checkboxGroupInput("ROITypes",label = "ROI Type",
                           choiceNames = list("TSS","CGI","Shores"),
                           choiceValues = list("TSS","CGI","SHO"),
                           selected = list("TSS","CGI")),
        
        actionButton("updateROIS","Update ROIS")),
        
        mainPanel(width = 9,
               plotOutput("grois",height="800px",width="1000px"))
        
      )}),
    
    tabPanel("Genes on Karyotype",
             plotOutput("gkaryo",width="800px",height="800px")),
  
    
    tabPanel("Somatic alterations",
             plotOutput("gsom",width="800px",height="800px")),
    
    tabPanel("Correlation graph A",
             plotOutput("gcorA")),
    
    tabPanel("Correlation graph B",
             plotOutput("gcorB"))
    
  ))
  


# SHINY SERVER -------


server <- function(input,output,session) {
  
  
  
  observeEvent(input$toGenes,{
    families <- strsplit(input$families,split="(\n|,)")[[1]] %>% gsub(" +","",.)
    familyNames <- paste(families,panther2[families,"name"],sep=": ")
    genes <- familyToGene(families)
    
    updateTextAreaInput(session,"genes",value=paste(genes,collapse="\n"))
    updateTextAreaInput(session,"familyNames",value=paste(familyNames,collapse="\n\n"))
    
  })
  
  observeEvent(input$clearFamilies,{
    updateTextAreaInput(session,"families",value="")
    updateTextAreaInput(session,"familyNames",value="")
  })
  
  observeEvent(input$toFamilies,{
    genelist <- strsplit(input$genes,split="(\n|,)")[[1]] %>% gsub(" +","",.)
    families <- geneToFamily(genelist)
    familyNames <- paste(families,panther2[families,"name"],sep=": ")
    
    updateTextAreaInput(session,"families",value=paste(families,collapse="\n"))
    updateTextAreaInput(session,"familyNames",value=paste(familyNames,collapse="\n\n"))
    
  })
  
  observeEvent(input$clearGenes,
               updateTextAreaInput(session,"genes",value=""))
  
  
  observeEvent(input$updateROIS, {
    
    # Select ROIS from genes
    
    families <- strsplit(input$families,split="(\n|,)")[[1]] %>% gsub(" +","",.)
    families <- families[families != ""]
    # families <- families[families %in% panther$hmmpanther] # remove families not listed in panther
    
    genes <- strsplit(input$genes,split="(\n|,)")[[1]] %>% gsub(" +","",.)
    # genes <- genes[genes %in% panther$hmmpanther] # remove genes not listed in panther
    
    ROITypes <- input$ROITypes
    
    rois <- subset(geneToProbe,hgnc_symbol %in% genes & ROIType %in% ROITypes)
    
    nmax <- sapply(list(families,genes,rois$ROI), function(x) length(unique(x))) %>% max
    
    output$grois <- renderPlot(
      height = 800+(nmax-20)*10,
      expr = roisGraph(families,genes,ROITypes) 
    )
    
     output$gkaryo <- renderPlot(
       expr = {
        kp <- plotKaryotype()
        rois2 <- rois[,c("hgnc_symbol","ROI","chr")] %>% unique
        rois2$center <- sapply(strsplit(rois2$ROI,split="[:-]"),function(x) as.numeric(x[3:4]) %>% mean)
        kpPlotMarkers(kp,chr=rois2$chr,x=rois2$center,y=.1,labels = rois2$hgnc_symbol,cex=1,srt=45)
       }
     )
    
  })
  
}


# Functions to be run on the Shiny Server


# Graph showing the families, genes and ROIs ----

roisGraph <- function(families,genes,ROITypes) {
  rois <- subset(geneToProbe,hgnc_symbol %in% genes & ROIType %in% ROITypes)
  rois <- rois[,c("hgnc_symbol","ensembl_transcript_id","ROI","chr")]
  rois <- unique(rois)
  rois$center <- sapply(strsplit(rois$ROI,split="(:|-)"),function(x) mean(as.numeric(x[3:4])))

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

# Somatic alterations -----






shinyApp(ui,server)
