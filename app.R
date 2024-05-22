# Install and load required packages if not already installed
if (!require("dipsaus")) install.packages("dipsaus")
if (!require("DT")) install.packages("DT")
if (!require("dplyr")) install.packages("dplyr")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("fgsea")) BiocManager::install("fgsea")
if (!require("fgsea")) install.packages("fgsea")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("ggpubr")) install.packages("ggpubr")
if (!require("gprofiler2")) install.packages("gprofiler2")
if (!require("emojifont")) install.packages("emojifont")
if (!require("plotly")) install.packages("plotly")
if (!require("shiny")) install.packages("shiny")
if (!require("shiny")) install.packages("shinydashboard")
if (!require("shinyjs")) install.packages("shinyjs")
if (!require("shinyWidgets")) install.packages("shinyWidgets")
if (!require("stringr")) install.packages("stringr")
if (!require("visNetwork")) install.packages("visNetwork")

# Load required packages
library(dipsaus)
library(DT)
library(dplyr)
library(fgsea)
library(ggplot2)
library(ggpubr)
library(gprofiler2)
library(emojifont)
library(plotly)
library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinyWidgets)
library(stringr)
library(visNetwork)

# Set directory to the src directory of the app so that the path of databases 
# is correctly resolved later (otherwise they won't be found)
setwd(getSrcDirectory(function(){})[1])

# Set the size limit of uploaded files
options(shiny.maxRequestSize = 50*1024^2)

# This css block is for the following reasons:
# - add top padding to rows to avoid overlapping of elements
# - move the length select box to the left to be more consistent
# - move the filter box to the left and add some bottom padding to be more consistent
# - move the filter box up a bit to align it with the length select box
css <- HTML("
/* customizing padding for buttons in the dashboard body */
  .row .pad-top {
     padding-top:25px;
  }
  
  /* customizing padding help button in the dashboard header */
  .main-header .dropdown {
     padding-top:9px;
     padding-right:25px;
  }
  
  /* customizing buttons of the DT table */
  .dataTables_wrapper .dataTables_length {
     float: left;
  }
  .dataTables_wrapper .dataTables_filter {
     float: left;
     padding-left: 50px;
     padding-bottom: 5px;
     margin-top: -2px;
  }
")

###################
#       UI       #
##################
ui <- fluidPage(

  # Set up css for the whole app
  tags$style(css),

  # Set up shinyjs so that we can use some of its functions
  # (e.g. hide, show, toggle)
  useShinyjs(), 
  
  # Use dashboard for the application
  dashboardPage(
    # HEADER
    dashboardHeader(title = span("EnrichRadar", style = "font-size: 40px"),
                    tags$li(class = "dropdown", actionButtonStyled(inputId = "HelpButton", 
                                                                   label = "Help", 
                                                                   type = "warning"))),
    
    # SIDEBAR
    dashboardSidebar(
      # This layout contains all the controls for filtering and downloading the data
      # Deseq2 output file input
      fileInput("csv", "Choose DESeq2 Output File", accept = c(".csv")),
      # Filtering by Log2FoldChange & padj
      sliderInput(inputId = "thresholdSliderLOG2FC",
                  label = "Log2FoldChange Threshold: ",
                  min = 0, max = 5, value = 0, step = 0.5),
      sliderTextInput(inputId = "thresholdSliderPADJ",
                      label = "padj Threshold: ",
                      choices = c(0.01, 0.05, "NONE"),
                      selected = "NONE",
                      grid = TRUE),
      # Filtering by DE type (down, up or both)
      radioButtons(inputId = "DEside",
                   label = "DE type:",
                   choices = c("Both" = "both",
                               "Down Regulated" = "down",
                               "Up Regulated" = "up"),
                   selected = "both"),
      # Filtering by coding genes or not
      checkboxInput(inputId = "KeepCodingGenes",
                    label = "Show only protein coding genes",
                    value = FALSE),
      # Filtering by known genes or not
      checkboxInput(inputId = "KeepKnownGenes",
                    label = "Show only known genes",
                    value = FALSE)
      ),
      # End of sidebar panel
    
    # BODY
    dashboardBody(
      ## Main Dataframe ##
      fluidRow(
        column(12, align = "center",
               h2(textOutput(outputId = "starterText")))
      ),
      # Title of the table (filename)
      h3(textOutput(outputId = "DFtitle")),
      # Save current selection button (to download as a csv file)
      downloadButton(outputId = "downloadMainTable",
                     label = "Save current table"),
      p(),
      # The main table
      DTOutput(outputId = "table"),
      hr(),
      
      fluidRow(
        column(12, align = "center",
               h3(textOutput(outputId = "wrongDATAmessage")))
      ),
      
      ## Volcano plot ##
      fluidRow(
        column(12, align = "center",
               h4(textOutput(outputId = "VolcanoMainText")))
      ),
      fluidRow(
        column(12, align = "center",
               actionButtonStyled(inputId = "launchVolcano", 
                                  label = "Volcano Plot", 
                                  type = "default")),
      ),
      p(),
      fluidRow(
        column(9, 
               uiOutput(outputId = "plotVolcano")),
        column(3, 
               plotOutput(outputId = "legendVolcano"))
      ),
      hr(),
      
      ## Gprofiler ##
      fluidRow(
        column(12, align = "center",
               h4(textOutput(outputId = "GprofilerMainText")))
      ),
      fluidRow(
        column(6, align = "right", class = "pad-top",
               actionButtonStyled(inputId = "launchGprofiler", 
                                  label = "Gprofiler", 
                                  type = "default")),
        column(6, align = "left",
               selectInput(inputId = "chooseOrganism", 
                           label = "Choose a specie",
                           choices = c("Human","Dog"),
                           selected = "Human",
                           width = "200px"))
      ),
      fluidRow(
        column(12, align = "center", 
               h5(textOutput(outputId = "GPtext"))),
      ),
      tabsetPanel(id = "TabsetGprofiler",
                  tabPanel(title = "Plot", 
                           fluidRow(
                             column(10, uiOutput(outputId = "plotPathways")),
                             column(2, selectInput(inputId = "nbTerm", 
                                                   tags$div("You have a lot of terms!", 
                                                            tags$br(),"Choose how many to display"),
                                                   choices = c("10","20","30","ALL"),
                                                   selected = "ALL",
                                                   width = "200px")))),
                  tabPanel(title = "Table", 
                           DTOutput(outputId = "gprofilerTable"),
                           downloadButton(outputId = "downloadGprofilerTable", 
                                          label = "Save Gprofiler table"))
      ),
      hr(),
      
      ## GSEA ##
      fluidRow(
        column(12, align = "center",
               h4(textOutput(outputId = "GSEAMainText"))),
      ),
      fluidRow(
        column(6, align = "right", class = "pad-top", 
               actionButtonStyled(inputId = "launchFGSEA", 
                                  label = "GSEA",
                                  type="default")),
        column(6, align = "left",
               selectInput(inputId = "chooseGSEADB",
                           label = "Choose a database",
                           choices = c("GO:MF","GO:CC","GO:BP","KEGG","REACTOME","WikiPathways","TRANSFAC & JASPAR PWMs","miRTarBase","CORUM","Human Phenotype Ontology"),
                           width = "200px")),
      ),
      fluidRow(
        column(12, align = "center", 
               h5(textOutput(outputId = "GSEAtext"))),
      ),
      tabsetPanel(id = "TabsetGSEA",
                  tabPanel(title = "Plot", 
                           plotOutput(outputId = "barplotGSEA")),
                  tabPanel(title = "Network",
                           p(), # Add an empty line for style
                           fluidRow(
                             column(12, align = "left",
                                    h5(textOutput(outputId = "GSEAtext2")))),
                           visNetworkOutput(outputId = "pathway_network")),
                  tabPanel(title = "Table", 
                           DTOutput(outputId = "fgseaTable"),
                           downloadButton(outputId = "downloadGSEATable",
                                          label = "Save GSEA table"))
      ),
      hr()
    )
  )
) 

###################
#     SERVER     #
##################
server <- function(input, output, session) {
  
  # Help button
  observeEvent(input$HelpButton, {
    
    showModal(modalDialog(
      # Volcano plot
      tags$h2('Volcano plot'),
      tags$h4('A volcano plot is a graphical tool used in gene expression analysis to visualize the relationship between fold change and statistical significance, 
              helping identify genes that are significantly differentially expressed between experimental conditions.'),
      # Gprofiler
      tags$h2('Gprofiler'),
      tags$h4('performs functional enrichment analysis, also known as over-representation analysis (ORA) or gene set enrichment analysis,
              on input gene list. It maps genes to known functional information sources and detects statistically significantly enriched terms'),
      # GSEA
      tags$h2('GSEA'),
      tags$h4('Gene Set Enrichment Analysis (GSEA) is a computational method that determines whether an a priori defined set of genes shows statistically 
      significant, concordant differences between two biological states (e.g. phenotypes). '),
      # Quit button
      footer = actionButton('QuitHelp', 'Ok'),
    ))
  })
  
  # When quit button is clicked on help, close the help window
  observeEvent(input$QuitHelp, {
    removeModal()
  })
  
  
  
  # Main text for Gprofiler analysis
  output$starterText <- renderText({
    "Click on \"Browse\" button of the side panel and choose your DESeq2 output file to start"
  })
  
  # Hide buttons until data is loaded
  shinyjs::hide("downloadMainTable")
  shinyjs::hide("launchVolcano")
  shinyjs::hide("launchGprofiler")
  shinyjs::hide("chooseOrganism")
  shinyjs::hide("launchFGSEA")
  shinyjs::hide("chooseGSEADB")
  
  # Hide tabset and volcano plot space until buttons are clicked
  shinyjs::hide("plotVolcano")
  shinyjs::hide("legendVolcano")
  shinyjs::hide("TabsetGprofiler")
  shinyjs::hide("TabsetGSEA")
  
  # Hide button to select how much term are displayed in gprofiler plot until condition is met (nbTerm > 50)
  shinyjs::hide("nbTerm")
  
  # Hide main/hint text
  shinyjs::hide("DFtitle")
  shinyjs::hide("wrongDATAmessage")
  shinyjs::hide("VolcanoMainText")
  shinyjs::hide("GprofilerMainText")
  shinyjs::hide("GSEAMainText")
  shinyjs::hide("GSEAtext")
  
  # Read CSV file
  data <- reactive({
    req(input$csv)
    read.csv(input$csv$datapath)
  })
  
  observe(
    if (!is.null(data())){
      # Hide starter text
      shinyjs::hide("starterText")
      # Show download/analysis buttons
      shinyjs::show("downloadMainTable")
      shinyjs::show("launchVolcano")
      shinyjs::show("launchGprofiler")
      shinyjs::show("chooseOrganism")
      shinyjs::show("launchFGSEA")
      shinyjs::show("chooseGSEADB")
      # Show main/hint text
      shinyjs::show("DFtitle")
      shinyjs::show("VolcanoMainText")
      shinyjs::show("GprofilerMainText")
      shinyjs::show("GSEAMainText")
      shinyjs::show("GSEAtext")
    }
  )
  
  # Create a reactive expression for filtered data
  filtered_data <- reactive({
    thresholdLOG2FC <- input$thresholdSliderLOG2FC
    thresholdPADJ <- input$thresholdSliderPADJ
    
    # Detect if the differential expression tool used is DESeq2 or edgeR
    if (all(c("log2FoldChange", "pvalue", "padj") %in% names(data()))){
      method = "deseq2"
    } else if (all(c("logFC", "PValue", "logCPM") %in% names(data()))){
      method = "edgeR"
    } else{
      # Not a differential expression file or not DESeq2 / edgeR
      method = "WRONG"
    }
    
    # Change the name of column to which the filters are applied depending if DESeq2 or edgeR
    if (method == "deseq2"){
      shinyjs::hide("wrongDATAmessage")
      # To filter dataframe with chosen thresholds in Log2FoldChange & padj with DESeq2 data
      if (input$DEside == "both") {
        filtered_data <- data()[abs(data()$log2FoldChange) >= thresholdLOG2FC & data()$padj <= thresholdPADJ, ]
      } else if (input$DEside == "up") {
        filtered_data <- data()[data()$log2FoldChange >= thresholdLOG2FC & data()$padj <= thresholdPADJ, ]
      } else if (input$DEside == "down") {
        filtered_data <- data()[data()$log2FoldChange <= 0-thresholdLOG2FC & data()$padj <= thresholdPADJ, ]
      }
    } else if (method == "edgeR"){
      shinyjs::hide("wrongDATAmessage")
      # To filter dataframe with chosen thresholds in Log2FoldChange & pvalue with edgeR data
      if (input$DEside == "both") {
        filtered_data <- data()[abs(data()$logFC) >= thresholdLOG2FC & data()$PValue <= thresholdPADJ, ]
      } else if (input$DEside == "up") {
        filtered_data <- data()[data()$logFC >= thresholdLOG2FC & data()$PValue <= thresholdPADJ, ]
      } else if (input$DEside == "down") {
        filtered_data <- data()[data()$logFC <= 0-thresholdLOG2FC & data()$PValue <= thresholdPADJ, ]
      }
      # rename first column as geneID if it is named X
      if (names(data()[1]) == "X"){
        filtered_data <- filtered_data %>% 
          rename(geneID = X)
      }
    } else if (method == "WRONG"){
      # If the data uploaded are not good (not deseq2 or edgeR, or something else entirely, print a error message and hide all the UI)
        filtered_data <- NULL
        output$wrongDATAmessage <- renderText({
          paste0(emoji("no_entry")," Currently, only files obtained with DESeq2 or edgeR are working with the application, please use a file obtained with one of those two tools and try again", emoji("no_entry"))
        })
        # Hide buttons
        shinyjs::hide("downloadMainTable")
        shinyjs::hide("launchVolcano")
        shinyjs::hide("launchGprofiler")
        shinyjs::hide("chooseOrganism")
        shinyjs::hide("launchFGSEA")
        shinyjs::hide("chooseGSEADB")
        # Hide main/hint text
        shinyjs::hide("DFtitle")
        shinyjs::hide("VolcanoMainText")
        shinyjs::hide("GprofilerMainText")
        shinyjs::hide("GSEAMainText")
        shinyjs::hide("GSEAtext")
        # Show data format error
        shinyjs::show("wrongDATAmessage")
        
        return(filtered_data)
    }
  
    
    # To avoid getting empty rows when moving thresholds: keeps only rows when the full line is not empty
    filtered_data <- filtered_data[rowSums(is.na(filtered_data)) != ncol(filtered_data), ]
    
    # Get the full dataframe even with non-expressed genes (with NA in Log2FoldChange and padj) when thresholds are by default
    if (thresholdLOG2FC == 0 & thresholdPADJ == "NONE") {
      # If edgeR was used replaced the X column by geneID
      if (names(data()[1]) == "X"){
        filtered_data <- data() %>%
          rename(geneID = X)
      } else {
        # If DESeq2 was used no need to replace column name
        filtered_data <- data()
      }
    }
    if (input$KeepCodingGenes == TRUE) {
      filtered_data <- filtered_data[filtered_data$gene_biotype == "protein_coding",]
    }
    if (input$KeepKnownGenes == TRUE) {
      filtered_data <- filtered_data[grepl("^(ENS|NM|NR)", filtered_data$geneID), ] # Keep only gene ID starting with prefix of Ensembl or RefSeq
    }
    
    return(filtered_data)
  })
  
  # Title to the table (file name)
  output$DFtitle <- renderText({
    paste0("File analyzed: ", input$csv$name)
  })
  
  # Render the filtered table
  output$table <- renderDT({
    datatable(filtered_data(), options = list(ordering = TRUE, pageLength = 10), rownames = FALSE)
  })
  
  # Download handler
  output$downloadMainTable <- downloadHandler(
    filename = function() {
      paste("DESeq2Viewer_LOG2FC_", input$thresholdSliderLOG2FC, "_PADJ_", input$thresholdSliderPADJ, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(filtered_data(), file, row.names = FALSE, quote = FALSE)
    }
  )
  
  
  
  # Main text for Volcano plot
  output$VolcanoMainText <- renderText({
    "Differential expression representation"
  })
  
  # Main text for Gprofiler analysis
  output$GprofilerMainText <- renderText({
    "Enrichment analysis with Gprofiler (analysis on several databases)"
  })
  
  # Main text for GSEA analysis
  output$GSEAMainText <- renderText({
    "Enrichment analysis with GSEA"
  })
  
  # Hint text for GSEA
  output$GSEAtext <- renderText({
    "NB: GSEA analysis is run on the whole data, you can't/there is no need to filter the data first"
  })
  
  
  
  ###################################
  # Volcano plot
  ###################################
  observeEvent(input$launchVolcano, {
    # Observe() so the plot changes when the table is updated
    observe({
      volcano_data <- select(data(), geneID,log2FoldChange,padj,gene_biotype,geneID,gene_name)
      volcano_data$padj <- ifelse((volcano_data$padj == 0 | -log(volcano_data$padj) > 30), exp(-30), volcano_data$padj)
      volcano_data$gene_biotype <- ifelse(is.na(volcano_data$gene_biotype), "other", volcano_data$gene_biotype)
      
      thresholdLOG2FC <- input$thresholdSliderLOG2FC
      thresholdPADJ <- input$thresholdSliderPADJ
      volcano_data$diff <- ifelse((volcano_data$log2FoldChange >= thresholdLOG2FC) & (volcano_data$padj <= thresholdPADJ),"UPPER", 
                                  ifelse((volcano_data$log2FoldChange <= 0-thresholdLOG2FC) & (volcano_data$padj <= thresholdPADJ), "UNDER", "NONE"))
      
      volcano_data <- filter(volcano_data, diff != "NONE")
      volcano_data <- volcano_data[complete.cases(volcano_data$padj,volcano_data$log2FoldChange,volcano_data$gene_biotype),]
      
      volcano_data$gene_annot <- ifelse(startsWith(volcano_data$geneID,"ENS"),"KNOWN", "NOVEL") 
      
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
        as_ggplot(VolcanoLegend)
      })
      
      shinyjs::show("plotVolcano")
      shinyjs::show("legendVolcano")
      
    })
  })
  
  
  ###################################
  # Gprofiler #
  ##################################
  observeEvent(input$launchGprofiler, {
    
    # Specie to use for gprofiler
    if (input$chooseOrganism == "Human"){
      specie <- "hsapiens"
    } else if (input$chooseOrganism == "Dog"){
      specie <- "clfamiliaris"
    }
    
    # Prepare and launch Gprofiler
    gprofilerdata <- filtered_data() %>%
      arrange(desc(abs(log2FoldChange)))
    
    GprofilerGeneList <- gprofilerdata[["geneID"]][1:nrow(gprofilerdata)] 
    
    # Run Gprofiler ONLY if the list of gene in query has less than 1000 genes, to avoid crash of the application
    if (length(GprofilerGeneList) <= 1000) {
      
      # Make the output of plot NULL so that it became blank and it shows that it is loading / doing something
      output$plotPathways <- renderUI({
        NULL
      })
        
      gostres <- gost(query = GprofilerGeneList, organism = specie, ordered_query = TRUE, evcodes = TRUE)
      
      # Reset pathways table and plot in case gostres returns nothing, to not keep the display of the last gprofiler run
      output$gprofilerTable <<- NULL
      output$plotPathways <<- NULL
      
      # Error message if bug
      if (is.null(gostres) == TRUE) {
        shinyjs::show("TabsetGprofiler")
        updateActionButtonStyled(session, "launchGprofiler", type="warning")
        output$plotPathways <- renderUI({
          HTML(
            as.character(div(style="color: orange", "No result to show, try with more genes / another set of genes"))
          )
        })   
        
        # Or make & print plot + table if there is no problem      
      } else {
        updateActionButtonStyled(session, "launchGprofiler", type="default")
        # Generate Gprofiler plot
        gost_plot <- gostplot(gostres, capped = FALSE, interactive = TRUE)
        
        # Table by pathways with all genes associated
        gostresDF <- as.data.frame(gostres$result)
        
        pathways_genes <- gostresDF %>% 
          dplyr::select(term_name, source, p_value, intersection)
        
        # Get number of genes associated with each term/pathway element
        gene_list <- pathways_genes$intersection
        gene_nb <- sapply(strsplit(gene_list,","), length)
        pathways_genes$gene_count <- gene_nb
        
        # Put a space after each coma in intersection column, to improve visualisation in dataframe
        pathways_genes$intersection <- gsub(",", ", ", pathways_genes$intersection)
        
        # Reorder dataframe columns
        pathways_genes <- pathways_genes %>%
          select(term_name, source, p_value, gene_count, intersection)
        
        # Reorder dataframe by pvalue (not descending order as the plot already reorder them)
        pathways_genes <- pathways_genes %>%
          arrange(p_value)
        
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
        
        # Download button
        output$downloadGprofilerTable <- downloadHandler(
          filename = function() {
            paste("Gprofiler_LOG2FC_", input$thresholdSliderLOG2FC, "_PADJ_", input$thresholdSliderPADJ, ".csv", sep = "")
          },
          content = function(file) {
            write.csv(pathways_genes, file, row.names = FALSE, quote = FALSE)
          }
        )
        
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
        if (nrow(pathways_genes) > 40) {
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
      
      # Make Gprofiler tabset appear 
      shinyjs::show("TabsetGprofiler")
      
      
      # Make Gprofiler tabset appear 
      output$GPtext <- renderText({
        paste0(emoji("warning"), " If you want to update the Gprofiler plot/table with a new selection of genes, you will need to press the button again. ", emoji("warning"))
      })
      
      shinyjs::show("GPtext")
      
    }
    
    # If there are more than 1000 genes in the query, display a warning message
    else {
      shinyjs::hide("nbTerm")
      shinyjs::show("TabsetGprofiler")
      updateActionButtonStyled(session, "launchGprofiler", type="warning")
      output$plotPathways <- renderUI({
        HTML(
          as.character(div(style="color: orange", "Query size too big. Max query size is 1000 genes. Please filter the genes before retrying."))
        )
      }) 
    }
    
  })

  
  
  ###################################
  # GSEA analysis with fgsea package
  ###################################

  observeEvent(input$launchFGSEA, {
    
    # Make the network blank to show that it is loading / doing something. Maybe also reduce lag
    output$pathway_network <- renderVisNetwork({
      NULL
    })
    
    # Log2FC for each genes (vector of Log2FC with genes ID associated as names)
    genes <- select(data(), gene_name, log2FoldChange)
    genes <- na.omit(genes)
    genes <- genes[order(genes$log2FoldChange, decreasing=TRUE),]
    
    ordered_gene_name <- genes$gene_name
    genes$gene_name <- NULL
    genes <- as.vector(genes$log2FoldChange)
    
    names(genes) <- ordered_gene_name
    

    # Import a pathway database and convert it into a table 
    processDB = function(DB){
      pathways=strsplit(scan(DB,sep="\n",what="character"),"\t")
      names(pathways)=sapply(pathways,function(x)x[1])
      pathways=lapply(pathways,function(x)x[3:length(x)])
      return(pathways)
    }
    
    # Changing the database used for GSEA analysis
    if (input$chooseGSEADB == "GO:MF"){
      file <- "GO_Molecular_Function_2023.txt"
    } else if (input$chooseGSEADB == "GO:CC"){
      file <- "GO_Cellular_Component_2023.txt"
    } else if (input$chooseGSEADB == "GO:BP"){
      file <- "GO_Biological_Process_2023.txt"
    } else if (input$chooseGSEADB == "KEGG"){
      file <- "KEGG_2021_Human.txt"
    } else if (input$chooseGSEADB == "REACTOME"){
      file <- "Reactome_2022.txt"
    } else if (input$chooseGSEADB == "WikiPathways"){
      file <- "WikiPathway_2023_Human.txt"
    } else if (input$chooseGSEADB == "TRANSFAC & JASPAR PWMs"){
      file <- "TRANSFAC_and_JASPAR_PWMs.txt"
    } else if (input$chooseGSEADB == "miRTarBase"){
      file <- "miRTarBase_2017.txt"
    } else if (input$chooseGSEADB == "CORUM"){
      file <- "CORUM.txt"
    } else if (input$chooseGSEADB == "Human Phenotype Ontology"){
      file <- "Human_Phenotype_Ontology.txt"
    }
    
    pathways <- processDB(paste0("./enrichR_databases/",file))
    
    
    # Define min and max number of genes associated with a pathway
    minSize = 3
    maxSize = 200
  
    # Launching GSEA analysis
    fgseaRes = data.frame(fgsea(pathways, genes, minSize=minSize, maxSize=maxSize))

    # Recovery of results
    cutpval = 0.05
    fdr = FALSE

    pv = ifelse(fdr == TRUE,"padj","pval")
    res = fgseaRes[which(fgseaRes[,pv] < cutpval),c("pathway", "padj","pval", "ES", "leadingEdge", "NES")]
    res$genes = unlist(lapply(res$leadingEdge, function(x) paste(x, collapse = ", ")))
    res = res[,c("pathway", pv, "ES","NES", "genes")]
    #print(head(res))
    respos = res[res[,"ES"]>0,]
    respos = respos[order(respos[,pv], decreasing = FALSE),]
    
    resneg = res[res[,"ES"]<0,]
    resneg = resneg[order(resneg[,pv], decreasing = FALSE),]
    
    # Barplot
    showCategory = 20
    respos_plot <- respos[1:min(showCategory, nrow(respos)),]
    resneg_plot <- resneg[1:min(showCategory, nrow(resneg)),]
    
    res_plot <- rbind(respos_plot, resneg_plot)
    res_plot[,"pval"]<- res_plot[,pv]
    
    # Make a column with shorter pathway names for pathway too long, to avoid excessive axis text size
    res_plot <- res_plot %>%
      mutate(pathway_shortname = ifelse(nchar(pathway) > 35, paste0(str_sub(pathway, end = 35), "..."), pathway))
    
    res_plot$pathway <- factor(res_plot$pathway, levels=res_plot$pathway)
    res_plot <- res_plot[order(res_plot[,"pval"], decreasing = TRUE),]
    
    barplot <- ggplot(res_plot, aes(x=reorder(pathway_shortname, -pval), y=-log10(pval), fill =NES))+
      geom_bar(stat='identity') + 
      scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red", limits=c(min(res_plot$NES), max(res_plot$NES))) +
      labs(x = "Pathway") +
      coord_flip() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12))
    
    
    # Pathway Network
    fgseaRes <- fgseaRes %>%
      filter(pval<0.05) # Keep only pathway significant (####change this later for res####)
    
    edges <- data.frame(from = character(), to = character(), NES = character())
    nodes <- data.frame(id = character(), group = character())
    
    for (line in 1:nrow(fgseaRes)){
      pathway <- fgseaRes[line,"pathway"]
      leadingEdge <- as.vector(fgseaRes[line,"leadingEdge"][[1]])
      for (gene in leadingEdge) {
        
        # add new edges
        new_row <- data.frame(from = gene, to = fgseaRes[line,"pathway"], NES = ifelse(fgseaRes[line,"NES"] > 0, "up", "down"))
        edges <- rbind(edges, new_row)
      }
    }
    
    # count number of edges for each genes
    edges_count <- edges %>%
      group_by(from) %>%
      mutate(edgesCount = n()) %>%
      select(-to, -NES)
    
    
    # add nodes + groups of nodes + number of edges for genes
    for (edge in unique(edges$from)){
      new_row <- data.frame(id = edge, group ="gene")
      nodes <- rbind(nodes, new_row)
    }
    
    for (edge in unique(edges$to)){
      new_row <- data.frame(id = edge, group ="pathway")
      nodes <- rbind(nodes, new_row)
    }
    
    nodes <- left_join(x = nodes, y = edges, join_by(id == to)) %>%
      select(-from) %>%
      unique()
    
    
    nodes <- left_join(x = nodes, y = edges_count, join_by(id == from)) %>%
      unique()
    
    nodes$group <- ifelse(nodes$group == "pathway", 
                          paste(nodes$group, nodes$NES, sep = "_"),
                          nodes$group)
    
    nodes$NES = NULL
    
    # add titles to nodes for caption
    nodes$title <- nodes$id
    
  
    # Network
    pathway_network <- visNetwork(nodes, edges) %>%
      visGroups(groupname = "gene", color = list(background = "palegreen", border = "seagreen", hover = list(background = "seagreen", border = "forestgreen"), highlight = list(background = "palegreen", border = "black")), size = 25) %>%
      visGroups(groupname = "pathway_up", color = list(background = "tomato", border = "firebrick", hover = list(background = "firebrick", border = "darkred"), highlight = list(background = "tomato", border = "black")), size = 50) %>%
      visGroups(groupname = "pathway_down", color = list(background = "lightblue", border = "royalblue", hover = list(background = "royalblue", border = "darkblue"), highlight = list(background = "lightblue", border = "black")), size = 50) %>%
      visLegend(addNodes = list(
        list(label = "Gene", shape = "circle", color = list(background = "palegreen", border = "seagreen")),
        list(label = "Pathway \n upregulated", shape = "circle", color = list(background = "tomato", border = "firebrick")),
        list(label = "Pathway \n downregulated", shape = "circle", color = list(background = "lightblue", border = "royalblue"))),
        useGroups = FALSE) %>%
      visOptions(highlightNearest = TRUE,
                 selectedBy = list(variable = "edgesCount", highlight = TRUE),
      ) %>%
      visInteraction(zoomView = TRUE, hover = TRUE, hoverConnectedEdges = TRUE, tooltipDelay = 0, navigationButtons = TRUE) %>%
      #visIgraphLayout() %>%
      visPhysics(stabilization = FALSE)

    
    
    # render GSEA Barplot
    output$barplotGSEA <- renderPlot({
      barplot
    })
    
    # render Pathway network
    output$pathway_network <- renderVisNetwork({
      pathway_network
    })
    
    # render the GSEA table
    output$fgseaTable <- renderDT({
      datatable(res, options = list(ordering = TRUE, order = list(1, 'asc'), pageLength = 5), rownames = FALSE)
    })
    
    # Download button
    output$downloadGSEATable <- downloadHandler(
      database_clean_name <- sub("\\.txt$", "", file),
      filename = function() {
        paste("GSEA_", database_clean_name, ".csv", sep = "")
      },
      content = function(file) {
        write.csv(res, file, row.names = FALSE, quote = FALSE)
      }
    )
    
    # Make GSEA tabset appear
    shinyjs::show("TabsetGSEA")
    
    # Second hint text for GSEA
    output$GSEAtext2 <- renderText({
      "Tips: Zoom in or hover on nodes to see the names of the genes/pathways"
    })
    
  })

}

# Run the application
shinyApp(ui, server)
