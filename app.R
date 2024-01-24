# Install and load required packages if not already installed
if (!require(shiny)) install.packages("shiny")
if (!require(DT)) install.packages("DT")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(plotly)) install.packages("plotly")
if (!require(gprofiler2)) install.packages("gprofiler2")
if (!require(dipsaus)) install.packages("dipsaus")

library(shiny)
library(DT)
library(ggplot2)
library(plotly)
library(gprofiler2)
library(dipsaus)

# Set the size limit of uploaded files
options(shiny.maxRequestSize = 50*1024^2)

# Define UI
ui <- fluidPage(
  
  titlePanel("DESeq2 Viewer"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("csv", "Choose DESeq2 Output File", accept = c(".csv")),
      sliderInput("thresholdSliderLOG2FC", "Log2FoldChange Threshold: ", min = 0, max = 5, value = 0, step = 0.5),
      sliderInput("thresholdSliderPADJ", "padj Threshold: ", min = 0.01, max = 1, value = 1, step = 0.01),
      radioButtons("DEside", "DE type:",
                   c("Both" = "both",
                     "Down" = "down",
                     "Up" = "up"
                     )),
      checkboxInput("KeepCodingGenes", "Show only protein coding genes", value = FALSE),
      checkboxInput("KeepKnownGenes", "Show only known genes", value = FALSE),
      downloadButton("downloadCSV", "Save current selection"),
      
      actionButtonStyled("launchGprofiler", "Gprofiler!", type="default"),
      
      width = 3
    ),
    
    mainPanel(
      DTOutput("table"),
      uiOutput("plot")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Read CSV file
  data <- reactive({
    req(input$csv)
    read.csv(input$csv$datapath)
  })
  
  # Create a reactive expression for filtered data
  filtered_data <- reactive({
    thresholdLOG2FC <- input$thresholdSliderLOG2FC
    thresholdPADJ <- input$thresholdSliderPADJ
    # To filter dataframe with chosen thresholds in Log2FoldChange & padj
    if(input$DEside == "both"){
      filtered_data <- data()[abs(data()$log2FoldChange) >= thresholdLOG2FC & data()$padj <= thresholdPADJ, ]
    }
    else if(input$DEside == "up"){
      filtered_data <- data()[data()$log2FoldChange >= thresholdLOG2FC & data()$padj <= thresholdPADJ, ]
    }
    else if(input$DEside == "down"){
      filtered_data <- data()[data()$log2FoldChange <= 0-thresholdLOG2FC & data()$padj <= thresholdPADJ, ]
    }
    
    # To avoid getting empty rows when moving thresholds: keeps only rows when the full line is not empty
    filtered_data <- filtered_data[rowSums(is.na(filtered_data)) != ncol(filtered_data), ]
    
    # Get the full dataframe even with non-expressed genes (with NA in Log2FoldChange and padj) when thresholds are by default
    if (thresholdLOG2FC == 0 & thresholdPADJ == 1){
      filtered_data <- data()
    }
    
    if (input$KeepCodingGenes == TRUE) {
      filtered_data <- filtered_data[filtered_data$gene_biotype == "protein_coding",]
    }
    
    if (input$KeepKnownGenes == TRUE) {
      filtered_data <- filtered_data[startsWith(filtered_data$geneID,"ENSG"),]
    }
    
    return(filtered_data)
  })
  
  # Render the filtered table
  output$table <- renderDT({
    datatable(filtered_data(), options = list(ordering = TRUE), rownames = FALSE)
  })
  
  # Download handler
  output$downloadCSV <- downloadHandler(
    filename = function() {
      paste("deseq2_output_LOG2FC_", input$thresholdSliderLOG2FC, "_PADJ_", input$thresholdSliderPADJ, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(filtered_data(), file, row.names = FALSE)
    }
  )
  
  # Gprofiler
  observeEvent(input$launchGprofiler, {
    GprofilerGeneList <- filtered_data()[["geneID"]][1:nrow(filtered_data())]
      
    gostres <- gost(query = GprofilerGeneList, organism = "hsapiens")
    
    output$plot <- renderUI({
        if (is.null(gostres) == TRUE){
          updateActionButtonStyled(session, "launchGprofiler", type="warning")
          HTML(
            as.character(div(style="color: orange", "No result to show, try with more genes / another set of genes"))
          )
        }
        else {
          gost_plot <- gostplot(gostres, capped = FALSE, interactive = TRUE)
          # Generate Gprofiler plot
          updateActionButtonStyled(session, "launchGprofiler", type="success")
          renderPlotly(
            gost_plot
          )
        }
    })
  })

  
}

# Run the application
shinyApp(ui, server)
