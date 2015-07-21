
# This is the server logic for a Shiny web application.
#
## To deploy, use:
#shinyapps::deployApp('~/Box Sync/habDX/habShiny/', appName ="habDX")

library(shiny)
library(dplyr)
library(ggplot2)


shinyServer(function(input, output) {
  
  geneNames <-read.table("MasterList_GeneNamesTable.txt", header=T, colClasses=c("character", "character"))
  geneNames<- mutate(geneNames, type= substr(Name,1,1))
  
  val <- reactiveValues()
  
  observeEvent(label = "newData",
               eventExpr = input$dataset,
               handlerExpr ={
                 withProgress(message = 'Loading data...',
                              value = 0, {for (i in 1:10) {incProgress(1/10)} })
                 val$fileName <- as.character(input$dataset);
                 val$currentData <- read.table(val$fileName, header=T)
               }
  )
  
  datReact <- reactive({   
    val$currentData %>%
      filter( !input$onlyITAG | geneNames$type[match(Gene, geneNames$Name)] == "S" ) %>% 
      filter(P.Value < 10^(-input$pval) ) %>% filter(input$includeDown | logFC > 0) %>% 
      mutate( logThres = factor(input$logFC < abs(logFC) )) 
  } )
  
  geneTable <-reactive({
    datReact()%>%
      filter(logThres == TRUE ) %>% 
      arrange(desc(abs(logFC) )) %>% 
      transmute(logFC = logFC, 
                AvgExpr = AveExpr, 
                p = P.Value, 
                Gene = Gene, 
                Note = geneNames$Note[match(Gene, geneNames$Name)])
  })
  
  output$histPlot <- renderPlot({
    ggplot(datReact()) + geom_bar(binwidth = 0.5, aes(x = logFC, fill = logThres )) +
      scale_fill_manual(values=c("grey90", "grey50"))+
      theme_classic()+theme(legend.position = "none")+
      geom_vline(xintercept = c(input$logFC, -input$logFC*as.numeric(input$includeDown)), colour= c("green", "red"))+
      labs(x="Expression difference A-B (log2 fold change)") 
  })
  
  output$selectedNumber <- renderText({
    selected <- sum( datReact()$logThres ==TRUE )
    total <- length(val$currentData[,1])
    sprintf("Selected %i genes (%.1f%% of total)", selected, 100*selected/total)
  })
  
  output$topTable <- renderTable({
    tt<-head(geneTable(), n=20)
    tt$p<-prettyNum(tt$p, format='e', digits=3)
    tt
  })
  
  output$downloadGeneTable <- downloadHandler(filename =
                                                function() {paste('Hab_genes_', Sys.Date(), '.txt', sep='')}, 
                                              content = 
                                                function(file) {write.table(geneTable(), file, row.names=F)})
  
  ####### Second tab
  normExp <- read.table("hab_flatNormalExpr.txt", header=T)
  treat <- read.table("hab_metadata.txt", header=T)

  observeEvent(eventExpr = input$singleGeneButton,
               label = "checkGene",
               handlerExpr ={
                 idx = match(input$userGene, row.names(normExp))
                 if(is.na(idx)){ val$buttonSwitch <- NaN
                                 val$errorMsg <- data.frame(Gene = input$userGene, 
                                                            Expressed ="No", 
                                                            Recognized_as_gene = !is.na(match(input$userGene, geneNames$Name)) )
                 } else { val$buttonSwitch <- 1
                          val$single <- mutate(treat, Acc=as.factor(Acc), expr=unlist(normExp[idx,]) )}
               })
  
  output$errorText <- renderTable({
    if(is.null(val$single)) return()
    if(is.na(val$buttonSwitch)) return(val$errorMsg)
  })
  
  output$singleGenePlot <- renderPlot({
    if(is.null(val$single)) return()
    if(is.na(val$buttonSwitch)) return()
    ggplot(val$single,aes(x=Tissue, y = expr, color=Acc ))+ scale_colour_discrete() + geom_point(size=3)
  })
}) #end shinyServer fn
