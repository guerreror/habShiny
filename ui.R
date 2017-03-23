
# This is the user-interface definition of a Shiny web application.

library(shiny)

dataset_names <- list(
  "Pollen (A) vs. Style (B)"= as.character("hab_pollen_vs_style_table.txt"),
  "Pollen IP (A) vs. no IP (B)" = as.character("hab_pollenIP_table.txt"),
  "Style UI (A) vs. no UI (B)" = as.character("hab_styleUI_table.txt"),
  "Style IP (A) vs. no IP (B)"= as.character("hab_styleIP_table.txt")
)

fluidPage(tabsetPanel(
  tabPanel("Lists of genes",
           pageWithSidebar(
             # Application title
             headerPanel("S.habrochaites differential gene expression"),
             
             # Sidebar with input
             sidebarPanel(
               # dataset drop-down
               selectInput("dataset", "Choose a contrast:", dataset_names),
               
               # pVal bar
               sliderInput("pval", "P-value threshold (-log10):", min = 2, max = 8, value = 3),
               
               # fold-change bar
               sliderInput("logFC", "Fold-change (log2):", min = 1, max = 15, value = 1),
               
               # checkbox for Lyco genes only
               checkboxInput("onlyITAG", "ITAG genes only (genes from S. lycopersicum)", TRUE),
               
               # include down-reg genes
               checkboxInput("includeDown", "Include down-regulated genes", FALSE),

               # download button
               downloadButton('downloadGeneTable', 'Download table of selected genes'),
               
               br(),
               br(),
               h4("Description of contrasts"),
               p(em("Pollen vs. Style:"),
                 "Genes different between tissues, across all hab accessions."),
               p(em("Pollen-side IP:"), 
                 "Genes different in pollen of accessions [1777,2119,1264] vs. [2098,407,1223], excludes genes that are also different in styles."),
               p(em("Style-side UI:"),
                 "Genes different in styles of accessions [1777,2098,2119] vs. [407,1223], excludes genes that are also different in pollen."),
               p(em("Style-side IP:"),
                 "Genes different in styles of accessions [1777, 2098] vs. [2119,1223], excludes genes that are also different in pollen.")
             ),
             
             #main for output
             mainPanel(
               plotOutput("histPlot"),
               h3(textOutput("selectedNumber")),
               p("Here is the Top-20:"),
               selectInput("sortBy", "Sort by ", list("Fold change" ="logFC", "Average expression"="AvgExpr")),
               tableOutput("topTable")
             )
           )
    ),
  
  tabPanel("Single gene trends",
           fluidPage(
             headerPanel("S.habrochaites differential gene expression"),
             
             fluidRow(
               column(width=4,
                      wellPanel(
                 textInput("userGene", label="Enter gene name", value="Solyc01g080970"),
                 actionButton("singleGeneButton", "Get gene expression plots")
                 )
               ),
               column(  width=7,
                 h3("Gene stats"),
                 p("From contrast selected in first tab. ", em("Rank") ," is the number of genes with logFC < focal gene,",em("percentile"),"is the fraction."),
                 tableOutput("singleRank")                 
                 )
                 ),
             fluidRow( 
               column(width=10,offset=1,
               tableOutput("errorText"),
               h3("Results across all samples"),
               plotOutput("singleGenePlot"),
               br(),
               h3("Pairwise comparisons"),
               br(),
               h4("Pollen-side IP plot:"),
               plotOutput("singleGene_poIP"),
               br(),
               h4("Style-side IP plot:"),
               plotOutput("singleGene_stIP"),
               br(),
               h4("Style-side UI plot:"),
               plotOutput("singleGene_UI")
               )
               )
           )
  ),
  
  tabPanel("Get a sequence",
            pageWithSidebar(
              headerPanel("S.habrochaites differential gene expression"),
              sidebarPanel(
                textInput("userSeq", label="Enter name (only 'c' genes)", value="c37227_g1_i1"),
                actionButton("singleSeqButton", "Get transcribed fragment sequence")
              ),
              mainPanel( 
                tableOutput("errorSeqChk"), 
                verbatimTextOutput("text.Seq")
                )
            )
 ),
 
 tabPanel("Search",
          fluidPage(
           headerPanel("S.habrochaites differential gene expression"),
           fluidRow(
             column(width=12,
             wellPanel(
             textInput("userSearch", label="Search in database:", value=""),
             checkboxInput("searchSelected", label="Include only genes selected in first tab", value=T ),
             actionButton("searchButton", "Go")
             )
             )
           ),
           fluidRow(
             tableOutput("searchResults")
             )
          )
 )
 
) #end Tabset Panel
) #end fluid page