# Install and load required packages if not already installed
if (!require(shiny)) install.packages("shiny")
if (!require(DT)) install.packages("DT")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(plotly)) install.packages("plotly")
if (!require(gprofiler2)) install.packages("gprofiler2")
if (!require(dipsaus)) install.packages("dipsaus")
if (!require(dplyr)) install.packages("dplyr")
if (!require(shinyjs)) install.packages("shinyjs")
if (!require(grid)) install.packages("grid")
if (!require(ggplotify)) install.packages("ggplotify")

library(shiny)
library(DT)
library(ggplot2)
library(plotly)
library(gprofiler2)
library(dipsaus)
library(dplyr)
library(shinyjs)
library(ggplotify)
library(grid)

# Set the size limit of uploaded files
options(shiny.maxRequestSize = 50*1024^2)

# Define UI
ui <- fluidPage(
  
  useShinyjs(),  # Set up shinyjs
  
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
      hr(),
      actionButtonStyled("launchGprofiler", "Gprofiler!", type="default"),
      actionButtonStyled("launchVolcano", "VolcanoPlot!", type="default"),
      
      width = 2
    ),
    
    mainPanel(
      DTOutput("table"),
      hr(),
      fluidRow(
        column(9, uiOutput("plotVolcano")),
        column(3, plotOutput("legendVolcano"))
      ),
      fluidRow(
        column(10, uiOutput("plotPathways")),
        column(2, selectInput("nbTerm", tags$div("You have a lot of terms!", tags$br(),"Choose how many to display"),
                              choices = c("20","30","40","ALL"),
                              selected = "ALL",
                              width = "200px")),
      ),
      DTOutput("gprofilerTable")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  shinyjs::hide("nbTerm")
  
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
    if (input$DEside == "both") {
      filtered_data <- data()[abs(data()$log2FoldChange) >= thresholdLOG2FC & data()$padj <= thresholdPADJ, ]
    } else if (input$DEside == "up") {
      filtered_data <- data()[data()$log2FoldChange >= thresholdLOG2FC & data()$padj <= thresholdPADJ, ]
    } else if (input$DEside == "down") {
      filtered_data <- data()[data()$log2FoldChange <= 0-thresholdLOG2FC & data()$padj <= thresholdPADJ, ]
    }
    
    # To avoid getting empty rows when moving thresholds: keeps only rows when the full line is not empty
    filtered_data <- filtered_data[rowSums(is.na(filtered_data)) != ncol(filtered_data), ]
    
    # Get the full dataframe even with non-expressed genes (with NA in Log2FoldChange and padj) when thresholds are by default
    if (thresholdLOG2FC == 0 & thresholdPADJ == 1) {
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
    datatable(filtered_data(), options = list(ordering = TRUE, pageLength = 10), rownames = FALSE)
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
      
    gostres <- gost(query = GprofilerGeneList, organism = "hsapiens", evcodes = TRUE)
    
    # Reset pathways table and plot in case gostres returns nothing, to not keep the display of the last gprofiler run
    output$gprofilerTable <<- NULL
    output$plotPathways <<- NULL
    
    # Plot
    output$plot <- renderUI({
        if (is.null(gostres) == TRUE) {
          updateActionButtonStyled(session, "launchGprofiler", type="warning")
          HTML(
            as.character(div(style="color: orange", "No result to show, try with more genes / another set of genes"))
          )
          
          
        } else {
          gost_plot <- gostplot(gostres, capped = FALSE, interactive = TRUE)
          # Generate Gprofiler plot
          updateActionButtonStyled(session, "launchGprofiler", type="success")
          #renderPlotly(
            #gost_plot
          #)
        }
    })
    
    if (!(is.null(gostres))) {
    
      # Table by pathways with all genes associated
      gostresDF <- as.data.frame(gostres$result)
      
      pathways_genes <- gostresDF %>% 
        dplyr::select(term_name, source, p_value, intersection)
      
      # Get number of genes associated with each term/pathway element
      gene_list <- pathways_genes$intersection
      gene_nb <- sapply(strsplit(gene_list,","), length)
      pathways_genes$gene_count <- gene_nb
      
      # Reorder dataframe
      pathways_genes <- pathways_genes %>%
        select(term_name, source, p_value, gene_count, intersection)
      
      # Print dataframe on screen
      output$gprofilerTable <- renderDT({
        datatable(pathways_genes, options = list(columnDefs = list(list(className = 'dt-center', targets = 3)),
                                                 order = list(2, 'asc'),
                                                 ordering = TRUE,
                                                 pageLength = 5),
                  rownames = FALSE) %>%
          formatStyle(
            target = 'row',
            columns = "source",
            backgroundColor = styleEqual(c("GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC", "TF", "MIRNA", "HPA", "CORUM", "HP", "WP"),
                                         c("#E361486D", "#46AB4B6D", "#FFAE486D", "#DD44776D", "#3366CC6D", "#5674A66D", "#23AB996D", "#6633CC6D", "#66AB016D", "#9A00996D", "#0099C66D")))
      })
    
    
      # Barplot with all the terms/pathways
      make_gprofiler_point <- function (pathways_genes_table) {
        BarplotPathways <- ggplot(pathways_genes_table, aes(x=p_value, y=reorder(term_name,-p_value), colour=source, size=gene_count, text=paste0("Term name: ",term_name, "<br>Source: ",source, "<br>p-value: ",p_value, "<br>Gene Count: ",gene_count))) +
          geom_point(stat="identity") +
          scale_colour_manual(name = "Source", values = c("GO:MF"="#E36148", "GO:CC"="#46AB4B", "GO:BP"="#FFAE48", "KEGG"="#DD4477", "REAC"="#3366CC", "TF"="#5674A6",
                                       "MIRNA"="#23AB99", "HPA"="#6633CC", "CORUM"="#66AB01", "HP"="#9A0099", "WP"="#0099C6")) +
          labs(x = "p-value", y = "", legend = "Source", title = "Enrichment terms/pathways", size = "", colour = "  Source") +
          theme(legend.text = element_text(size = 9),
                legend.title = element_text(size = 14),
                axis.title.x = element_text(size = 14),
                axis.text = element_text(size = 9) 
                )
      
        plotlyPlot <- ggplotly(BarplotPathways + scale_x_reverse(),
                          tooltip = "text") # Control what is display in the point label
        
        plotlyPlot <- plotlyPlot %>%
          layout(legend = list(itemsizing = "constant")) # Control size of points in legend
                 
        return (plotlyPlot)
      }
      
     
      # Button to reduce number of genes on the enrichment terms/pathways barplot if there are too many
      if (nrow(pathways_genes) > 50) {
        shinyjs::show("nbTerm")
        observe ({
          output$plotPathways <- renderUI({ 
            renderPlotly ({
              nbTerm_value <- as.integer(input$nbTerm)
              if (input$nbTerm == "ALL") {
                nbTerm_value = nrow(pathways_genes)
              }
              make_gprofiler_point(pathways_genes[1:nbTerm_value,])
            })
          })
        })
      } else {
        shinyjs::hide("nbTerm")
        output$plotPathways <- renderUI({
          renderPlotly ({
            make_gprofiler_point(pathways_genes)
          })
        })
      }
    }
  })
  
  
  # Volcano plot
  observeEvent(input$launchVolcano, {
    
    volcano_data <- select(data(), log2FoldChange,padj,gene_annot,gene_biotype,geneID,gene_name)
    volcano_data$padj <- ifelse((volcano_data$padj == 0 | -log(volcano_data$padj) > 30), exp(-30), volcano_data$padj)
    volcano_data$gene_biotype <- ifelse(is.na(volcano_data$gene_biotype), "other", volcano_data$gene_biotype)
    
    thresholdLOG2FC <- input$thresholdSliderLOG2FC
    thresholdPADJ <- input$thresholdSliderPADJ
    volcano_data$diff <- ifelse((volcano_data$log2FoldChange >= thresholdLOG2FC) & (volcano_data$padj <= thresholdPADJ),"UPPER", 
                                ifelse((volcano_data$log2FoldChange <= 0-thresholdLOG2FC) & (volcano_data$padj <= thresholdPADJ), "UNDER", "NONE"))
    
    volcano_data <- filter(volcano_data, diff != "NONE")
    volcano_data <- volcano_data[complete.cases(volcano_data$padj,volcano_data$log2FoldChange,volcano_data$gene_annot,volcano_data$gene_biotype),]
    
    VolcanoPlot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log(padj), shape = gene_biotype, fill = factor(diff), size = gene_annot, 
                                            text = paste0("GeneID: ",geneID, "<br>Gene name: ",gene_name, "<br>Log2FoldChange: ",log2FoldChange, "<br>p-adjusted: ",padj))) +
    geom_point(aes(stroke = .2)) +
    # add some lines
    geom_vline(xintercept = thresholdLOG2FC, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 0-thresholdLOG2FC, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 0, color = "black") +
    labs(x = "Log2FC", y = "-Log(p.adj)", fill="Differential Expression", size="Origin", shape="Gene Biotype") +
    guides(fill = guide_legend(override.aes = list(size = 6, shape=21)), shape = guide_legend(override.aes = list(size = 6))) +
    scale_fill_manual(values = c("UNDER"="#56B4E9", "UPPER"="#D55E00")) +
    scale_shape_manual(values = c("lncRNA" = 21 ,"protein_coding" = 23,"other" = 22)) +
    theme(legend.text = element_text(size = 16),
          legend.title = element_text(face = "bold", size = 18))
    
    # Extract legend to keep the ggplot type of legend
    VolcanoLegend <- get_legend(VolcanoPlot) 
    
    # Declare limits of the plot & delete legend 
    VolcanoPlot <- VolcanoPlot + scale_x_continuous(limits = c(-max(abs(volcano_data$log2FoldChange)), max(abs(volcano_data$log2FoldChange)))) +
      theme(legend.position="none")
    
    # Transform into plotly object
    VolcanoPlot <- ggplotly(VolcanoPlot,
                            tooltip = "text")
    
    # Print plot
    output$plotVolcano <- renderUI({
      renderPlotly({VolcanoPlot})
    })
    # Print legend
    output$legendVolcano <- renderPlot({
      grid.draw(VolcanoLegend)
    })
    
  })

}
# Run the application
shinyApp(ui, server)
