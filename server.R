
# This is the server logic for a Shiny web application.
#
## To deploy, use:
#shinyapps::deployApp('~/Box Sync/habDX/habShiny/', appName ="habDX")

library(shiny)
library(dplyr)
library(ggplot2)
library(readr)

shinyServer(function(input, output) {
  
  geneNames <-read_delim("MasterList_GeneNamesTable.txt",col_names = T, delim = "\t" )
  geneNames<- mutate(geneNames, Name= sapply(Name, function(x) strsplit(x, fixed=T, split=".")[[1]][1]),
                     Note=Note,
                     type= substr(Name,1,1))
  
  val <- reactiveValues()
  
  observeEvent(label = "newData",
               eventExpr = input$dataset,
               handlerExpr ={
                 withProgress(message = 'Loading data...',
                              value = 0, {for (i in 1:10) {incProgress(1/10)} })
                 val$fileName <- as.character(input$dataset);
                 val$currentData <- read_delim(val$fileName,col_names = T, delim =" ")
                 val$currentData$Gene <- sapply(val$currentData$Gene, function(x) strsplit(x, fixed=T, split=".")[[1]][1])
               }
  )
  
  datReact <- reactive({  
    withProgress(message = 'Loading data...',
                 value = 0, {for (i in 1:10) {incProgress(1/10)} })
    val$currentData %>%
      filter( !input$onlyITAG | geneNames$type[match(Gene, geneNames$Name )] == "S" ) %>% 
      filter(P.Value < 10^(-input$pval) ) %>% filter(input$includeDown | logFC > 0) %>% 
      mutate( logThres = factor(input$logFC < abs(logFC) )) 
  } )
  
  geneTable <-reactive({
    withProgress(message = 'Loading data...',
                 value = 0, {for (i in 1:10) {incProgress(1/10)} })
    out<-datReact()%>%
      filter(logThres == TRUE ) %>% 
      transmute(logFC = logFC, 
                AvgExpr = AveExpr, 
                p = P.Value, 
                Gene = Gene, 
                Note = geneNames$Note[match(Gene, geneNames$Name)])
    arrange(out, desc(get(input$sortBy, out) ))     
  })
  
  output$histPlot <- renderPlot({
    withProgress(message = 'Loading plot...',
                 value = 0, {for (i in 1:10) {incProgress(1/10)} })
    
    ggplot(datReact()) + geom_bar(binwidth = 0.5, aes(x = logFC, fill = logThres )) +
      scale_fill_manual(values=c("grey90", "grey50"))+
      theme_classic()+theme(legend.position = "none")+
      geom_vline(xintercept = c(input$logFC, -input$logFC*as.numeric(input$includeDown)), colour= c("green", "red"))+
      labs(x="Expression difference A-B (log2 fold change)") 
  })
  
  output$selectedNumber <- renderText({
    selected <- sum( datReact()$logThres ==TRUE )
    total <- dim(val$currentData)[1]
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
  normExp <- read_delim("hab_flatNormalExpr.txt", col_names = T, delim = " ")
  normExp$Gene <- sapply(normExp$Gene, function(x) strsplit(x, fixed=T, split=".")[[1]][1])
  treat <- read_delim("hab_metadata.txt", col_names=T, delim = " " )
  
  observeEvent(eventExpr = input$singleGeneButton,
               label = "checkGene",
               handlerExpr ={
                 withProgress(message = 'Loading data...',
                              value = 0, {for (i in 1:10) {incProgress(1/10)} })
                 
                 idx = match(input$userGene, normExp$Gene)
                 if(is.na(idx)){ val$buttonSwitch <- NaN
                                 val$errorMsg <- data.frame(Gene = input$userGene, 
                                                            Expressed ="No", 
                                                            Recognized_as_gene = !is.na(match(input$userGene, geneNames$Name)) )
                 } else { val$buttonSwitch <- 1
                          val$single <- mutate(treat, Acc=as.factor(Acc), expr=unlist(normExp[idx,-1]) )}
               })
  
  output$errorText <- renderTable({
    if(is.null(val$single)) return()
    if(is.na(val$buttonSwitch)) return(val$errorMsg)
  })
  
  output$singleGenePlot <- renderPlot({
    if(is.null(val$single)) return()
    if(is.na(val$buttonSwitch)) return()
    ggplot(val$single,aes(x=Tissue, color=Tissue, y = expr ))+
      facet_wrap(~Acc, nrow=1, scales="free_x") + 
      geom_point(size=3) + 
      theme_classic() + 
      theme(text=element_text(size=18), legend.position="none", panel.border=element_rect(fill=NA, color="black"), plot.background = element_rect(fill="grey95"))+
      labs(x = "", y="Normalized expression (log2)")  
  })
  
  output$singleRank <- renderTable({
    if(is.null(val$single)) return()
    if(is.na(val$buttonSwitch)) return()
    outS <- val$currentData[match(input$userGene, val$currentData$Gene),]
    out <- cbind(outS[,c(1,2,4)], sum(val$currentData$logFC <= outS$logFC), mean(val$currentData$logFC <= outS$logFC))
    names(out)[c(4,5)] <- c("Rank", "Percentile")
    out
  }, digits = 3, include.rownames = FALSE )
  
  output$singleGene_poIP <- renderPlot({
    if(is.null(val$single)) return()
    if(is.na(val$buttonSwitch)) return()
    poIPsingle <- filter(val$single, Tissue =="pollen" & poIP !="out")
    ggplot(poIPsingle)+ 
      geom_boxplot(aes(x = poIP, y = expr), fill="grey90") + 
      theme_classic()+ 
      theme(text=element_text(size=18),panel.border=element_rect(fill=NA, color="black", size=1),plot.background = element_rect(fill="grey95"))+ 
      labs(x = "pollen-side IP factor", y="Normalized expression (log2)")
  })
  
  output$singleGene_stIP <- renderPlot({
    if(is.null(val$single)) return()
    if(is.na(val$buttonSwitch)) return()
    stIPsingle <- filter(val$single, Tissue =="style" & stIP !="out")
    ggplot(stIPsingle)+ 
      geom_boxplot(aes(x = stIP, y = expr), fill="grey90") +
      theme_classic()+ theme(text=element_text(size=18),panel.border=element_rect(fill=NA, color="black", size=1),plot.background = element_rect(fill="grey95"))+ 
      labs(x = "style-side IP factor", y="Normalized expression (log2)")
  })
  
  output$singleGene_UI <- renderPlot({
    if(is.null(val$single)) return()
    if(is.na(val$buttonSwitch)) return()
    UIsingle <- filter(val$single, Tissue =="style" & UI !="out")
    ggplot(UIsingle)+ 
      geom_boxplot(aes(x = UI, y = expr), fill="grey90") +
      theme_classic()+ theme(text=element_text(size=18),panel.border=element_rect(fill=NA, color="black", size=1),plot.background = element_rect(fill="grey95"))+ 
      labs(x = "style-side UI factor", y="Normalized expression (log2)")
  })
  
  
  ###### Third tab
  seqs <- read_delim( "hE_transfrag_seqs_singleline.txt", col_names=F, delim = "\t" )
  names(seqs) <- c("gene", "sequence")
  observeEvent(eventExpr = input$singleSeqButton,
               label = "checkSeq",
               handlerExpr ={
                 idx = match(input$userSeq, seqs$gene)
                 if(is.na(idx)){ val$seq <- NaN
                                 val$seq_errorMsg <- data.frame(Gene = input$userSeq, 
                                                                Recognized_as_C_gene ="No", 
                                                                Recognized_as_gene = !is.na(match(input$userSeq, geneNames$Name)) )
                 } else { val$seq <- seqs$sequence[idx]}
               })
  
  output$errorSeqChk <- renderTable({
    if(is.null(val$seq)) return()
    if(is.na(val$seq)) return(val$seq_errorMsg)
  })
  
  output$text.Seq <- renderText({
    if(is.null(val$seq)) return()
    if(is.na(val$seq)) return("No sequence found.")
    val$seq
  })
  
  
  ##### Fourth tab
  
  searchTable <- eventReactive(
    input$searchButton,{
      withProgress(message = 'Searching...',
                   value = 0, {for (i in 1:10) {incProgress(1/10)} })
      if(!input$searchSelected){
        tab <-geneNames; cols <- c(1,2)
      }else{ tab <- geneTable()%>%mutate(Name=Gene, Gene=NULL);cols <- c(5,4,1,2,3) }
      
      got <- grepl(input$userSearch, tab$Name) | grepl(input$userSearch, tab$Note)
      out <- tab[got,cols]
      out
      
      })
  
  output$searchResults <- renderTable( {
    withProgress(message = 'Searching...',
                 value = 0, {for (i in 1:10) {incProgress(1/10)} })
    searchTable()
    }, include.rownames = FALSE)

}) #end shinyServer fn
