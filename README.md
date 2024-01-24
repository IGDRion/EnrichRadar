# DESeq2 viewer

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Description
This R shiny application allows to view the output of DESeq2 and make interactive filters over the data in an easy way. It is also able to produce an enrichment analysis via Gprofiler2 and return an interactive figure.

## Installation
If you have it, simply open the script in R studio and click on "Run App" on the top right corner of the source panel. You can also run it using the command **runApp('path/to/app.R')**. If you don't have R studio, you can open a R session in a terminal, install & load or load the library "shiny", and then run the command **runApp('path/to/app.R')**. It should open the application on your browser.

## Usage
To load a DESeq2 output file, click on the browse button and select the file you would like to view.
Once the file is loaded, you can filter your data :
- Choose wich threshold you want to use for Log2FoldChange and p-adujsted to filter your data
- Select which type of differential expression you want to see : "Both" shows down and up regulated genes, "Down" shows only down regulated genes, and "Up" shows only up regulated genes
- Two checkboxes allows you to filter if you want to see only protein coding genes & only known genes
- "Save current selection" button allows to save the dataframe with all the modification done to it (sorting of rows are not saved however)
- Gprofiler button generate an interactive graph with gprofiler2 (tool for multi-databases enrichment analysis) 

## Support
If you are having any trouble with the application, contact **victor.lebars@univ-rennes1.fr**

## Authors and acknowledgment
Author: Victor Le Bars

## License
This project is licensed under the [MIT License](LICENSE) - see the [LICENSE](LICENSE) file for details.