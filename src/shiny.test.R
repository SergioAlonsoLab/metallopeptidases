library(shiny)
library(dplyr)

setwd("/imppc/labs/mplab/share/metallopeptidases/")
load("data/pantherdb.rda",verbose = T)


ui <- fluidPage(
  tabsetPanel(
    tabPanel("Genes and Families", {
      fluidRow(
        
        column(width = 3,
               textAreaInput("genes","Genes",height = "500px",resize="none"),
               column(width = 6, actionButton("toFamilies","To Families")),
               column(width = 6, actionButton("clearGenes","Clear"))
        ),
        
        column(width = 3,
               textAreaInput("families","Families",height = "500px",resize="none"),
               column(width = 6,actionButton("toGenes","To Genes")),
               column(width = 6,actionButton("clearFamilies","Clear"))
        ),
        
        column(width = 6,
               textAreaInput("familyNames","Family Names",height = "500px",resize="none"))
      )
    }),
    
    tabPanel("ROIS", {
      fluidRow(
        column(width = 3,
               checkboxGroupInput("ROItypes",label = "ROI Type",
                           choiceNames = list("TSS","CGI","Shores"),
                           choiceValues = list("TSS","CGI","SHO"),
                           selected = list("TSS","CGI","SHO")),
        
        actionButton("updateROIS","Update ROIS"))
        
      )}
    )
  ))
  
  
server <- function(input,output,session) {
  
  
  
  observeEvent(input$toGenes,{
    families <- strsplit(input$families,split="(\n|,)")[[1]]
    familyNames <- paste(families,panther2[families,"name"],sep=": ")
    genes <- familyToGene(families)
    
    updateTextAreaInput(session,"genes",value=paste(genes,collapse="\n"))
    updateTextAreaInput(session,"familyNames",value=paste(familyNames,collapse="\n"))
    
  })
  
  observeEvent(input$clearFamilies,{
    updateTextAreaInput(session,"families",value="")
    updateTextAreaInput(session,"familyNames",value="")
  })
  
  observeEvent(input$toFamilies,{
    genelist <- strsplit(input$genes,split="(\n|,)")[[1]]
    families <- geneToFamily(genelist)
    familyNames <- paste(families,panther2[families,"name"],sep=": ")
    
    updateTextAreaInput(session,"families",value=paste(families,collapse="\n"))
    updateTextAreaInput(session,"familyNames",value=paste(familyNames,collapse="\n"))
    
  })
  
  observeEvent(input$clearGenes,
               updateTextAreaInput(session,"genes",value=""))
  
}

shinyApp(ui,server)
