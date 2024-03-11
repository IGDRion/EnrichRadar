# Install and load required packages if not already installed
if (!require(shiny)) install.packages("shiny")
if (!require(DT)) install.packages("DT")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(plotly)) install.packages("plotly")
if (!require(gprofiler2)) install.packages("gprofiler2")
if (!require(dipsaus)) install.packages("dipsaus")
if (!require(dplyr)) install.packages("dplyr")
if (!require(shinyjs)) install.packages("shinyjs")
if (!require(ggpubr)) install.packages("ggpubr")
if (!require(fgsea)) install.packages("fgsea")
if (!require(stringr)) install.packages("stringr")
if (!require(shinyWidgets)) install.packages("shinyWidgets")

library(shiny)
library(DT)
library(ggplot2)
library(plotly)
library(gprofiler2)
library(dipsaus)
library(dplyr)
library(shinyjs)
library(ggpubr)
library(fgsea)
library(stringr)
library(shinyWidgets)

# Set the size limit of uploaded files
options(shiny.maxRequestSize = 50*1024^2)

css <- HTML("
  .row .pad-top {
     padding-top:25px;
  }
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

# Define UI
ui <- fluidPage(
  tags$style(css), # Set up css
  
  useShinyjs(),  # Set up shinyjs
  
  titlePanel("DESeq2 Viewer"),

  
  sidebarLayout(
    sidebarPanel(
      fileInput("csv", "Choose DESeq2 Output File", accept = c(".csv")),
      sliderInput("thresholdSliderLOG2FC", "Log2FoldChange Threshold: ", min = 0, max = 5, value = 0, step = 0.5),
      sliderTextInput("thresholdSliderPADJ", "padj Threshold: ", choices = c(0.01, 0.05, "NONE"), selected = "NONE", grid = TRUE),
      radioButtons("DEside", "DE type:",
                   c("Both" = "both",
                     "Down" = "down",
                     "Up" = "up"
                     )),
      checkboxInput("KeepCodingGenes", "Show only protein coding genes", value = FALSE),
      checkboxInput("KeepKnownGenes", "Show only known genes", value = FALSE),
      downloadButton("downloadMainTable", "Save current selection"),
      hr(),
      
      width = 2
    ),
    
    mainPanel(
      DTOutput("table"),
      hr(),
      fluidRow(
        column(12, actionButtonStyled("launchVolcano", "Volcano Plot", type="default"), align = "center"),
      ),
      fluidRow(
        column(9, uiOutput("plotVolcano")),
        column(3, plotOutput("legendVolcano"))
      ),
      hr(),
      fluidRow(
        column(6, actionButtonStyled("launchGprofiler", "Gprofiler", type="default"), align = "right", class = "pad-top"),
        column(6,selectInput("chooseOrganism", "Please choose a specie",
                             choices = c("Human","Dog"),
                             selected = "Human",
                             width = "200px"), align = "left"),
      ),
      tabsetPanel(id="TabsetGprofiler",
        tabPanel("Plot", fluidRow(
          column(10, uiOutput("plotPathways")),
          column(2, selectInput("nbTerm", tags$div("You have a lot of terms!", tags$br(),"Choose how many to display"),
                                choices = c("20","30","40","ALL"),
                                selected = "ALL",
                                width = "200px")))),
        tabPanel("Table", DTOutput("gprofilerTable"),
                 downloadButton("downloadGprofilerTable", "Save Gprofiler table"))
      ),
      hr(),
      fluidRow(
        column(6,actionButtonStyled("launchFGSEA", "GSEA", type="default"), align = "right", class = "pad-top"),
        column(6,selectInput("chooseGSEADB", "Please choose a database",
                             choices = c("GO:MF","GO:CC","GO:BP","KEGG","REACTOME","WikiPathways","TRANSFAC & JASPAR PWMs","miRTarBase","CORUM","Human Phenotype Ontology"),
                             width = "200px"), align = "left"),
      ),
      tabsetPanel(id="TabsetGSEA",
        tabPanel("Plots", fluidRow(
          column(6, plotOutput("barplotGSEA")),
          column(6, plotOutput("lineplotGSEA"))
        )),
        tabPanel("Table", DTOutput("fgseaTable"),
                 downloadButton("downloadGSEATable", "Save GSEA table"))
      ),
      
      hr()
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Hide buttons until data is loaded
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
  
  # Read CSV file
  data <- reactive({
    req(input$csv)
    read.csv(input$csv$datapath)
  })
  
  observe(
    if (!is.null(data())){
      shinyjs::show("launchVolcano")
      shinyjs::show("launchGprofiler")
      shinyjs::show("chooseOrganism")
      shinyjs::show("launchFGSEA")
      shinyjs::show("chooseGSEADB")
    }
  )
  
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
    if (thresholdLOG2FC == 0 & thresholdPADJ == "NONE") {
      filtered_data <- data()
    }
    if (input$KeepCodingGenes == TRUE) {
      filtered_data <- filtered_data[filtered_data$gene_biotype == "protein_coding",]
    }
    if (input$KeepKnownGenes == TRUE) {
      filtered_data <- filtered_data[startsWith(filtered_data$geneID,"ENS"),]
    }
    
    return(filtered_data)
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
  
  
  
  ###################################
  # Volcano plot
  ###################################
  observeEvent(input$launchVolcano, {
    
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
    GprofilerGeneList <- filtered_data()[["geneID"]][1:nrow(filtered_data())]
      
    gostres <- gost(query = GprofilerGeneList, organism = specie, evcodes = TRUE)
    
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
    
    # Make Gprofiler tabset appear 
      shinyjs::show("TabsetGprofiler")
  })

  
  
  ###################################
  # GSEA analysis with fgsea package
  ###################################
  observeEvent(input$launchFGSEA, {
    
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
    
    pathways <- processDB(paste0("/Users/Victor/",file))
    
    
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
    
    p <- ggplot(res_plot, aes(x=reorder(pathway_shortname, -pval), y=-log10(pval), fill =NES))+
      geom_bar(stat='identity') + 
      scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red", limits=c(min(res_plot$NES), max(res_plot$NES))) +
      labs(x = "Pathway") +
      coord_flip() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12))
    
    
    # GSEA plot
    pathway_selection <- as.character(res_plot[order(res_plot$pval, decreasing = FALSE),"pathway"])
    p2 <- plotGseaTable(pathways[pathway_selection], genes, fgseaRes, gseaParam=0.5,
                        #axisLabelStyle = list(size=12),
                        pathwayLabelStyle = list(size=12))
    
    # render GSEA Barplot
    output$barplotGSEA <- renderPlot({
      p
    })
    
    # render GSEA plot
    output$lineplotGSEA <- renderPlot({
      p2
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
    
    # Make Gprofiler tabset appear
    shinyjs::show("TabsetGSEA")
    
  })

}
# Run the application
shinyApp(ui, server)
