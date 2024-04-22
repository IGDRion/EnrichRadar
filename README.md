# DESeq2 viewer

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Description
This R shiny application allows to view the output of DESeq2 and make interactive filters over the data in an easy way. It is also able to produce an enrichment analysis via Gprofiler2 and/or GSEA, and return interactive plots.

## Installation

- Git clone the repository:

```bash
git clone https://github.com/IGDRion/deseq2-viewer.git
```

- Install the required R libraries using e.g. with `conda`

```bash
conda install conda-forge::r-base conda-forge::r-shiny r::r-dt conda-forge::r-ggplot2  plotly::plotly conda-forge::r-gprofiler2 r-dipsaus conda-forge::r-dplyr conda-forge::r-shinyjs conda-forge::r-ggpubr bioconda::bioconductor-fgsea conda-forge::r-stringr
```

- Launch the app:

```bash
Rscript /path/to/app.R
```

It should open the application on your web browser : `http://127.0.0.1:7449`.

## Usage
### Loading Data
To load a DESeq2 output file, click on the `Browse...` button on the top left corner and select the file you would like to analyse. 

### Data Filtering
Once the file is loaded, the dataframe should now appear on the main page. With the panel on the left side, you can filter the table by several criteria:

- You can use the `Log2FoldChange Threshold` slider if you want to filter genes based on their differential expression level (By default all genes are kept).
- In the same way, you can use the `padj Threshold` slider to filter genes that have a padj value above the chosen threshold (By default all genes are kept).
- `DE type` buttons are to select which type of differential expression you want to see : "Both" shows down and up regulated genes (positive & negative Log2FoldChange), "Down" shows only down regulated genes (negative Log2FoldChange), and "Up" shows only up regulated genes (positive Log2FoldChange).
- Two checkboxes allows you to filter if you want to see only protein coding genes & only known genes.
- `Save current selection` button allows to save the dataframe with all the filters applied to it (sorting of rows are not saved however).

### Volcano plot / Enrichment analysis
After loading the file, three buttons appear on the main page:

- `Volcano Plot` generates a volcano plot. The data points displayed correspond with the filters applied on the side panel, like for the main table. You can see information about each point by hovering your mouse on it. It should update whenever you change your filters on the table.
- `Gprofiler` launch a query to several pathways databases with the current selection of gene. It generate a plot and a pable, accessible with two tabs. You can also select the organism on which you would like to submit the query (for now only Human and Dog are implemented in the application). Unlike the Volcano plot, you should click the button if you want to launch gprofiler again.
- `GSEA` launch a query to one specified pathway database, which you can select next to the button. Unlike the volcano plot or the Gprofiler analysis, it always launch on the complete list of genes, and it takes into account the Log2FoldChange of every gene into the analysis. Three tabs are available for this analysis: one barplot with over- and under- expressed pathways, one with a network representation of these pathways and genes associated (this can lag if too many nodes), and one with a table.


## Example file
A Toy example is available to try the application : `ToyExample.csv`.

The data of this example file are from the [**airway** dataset](https://bioconductor.org/packages/release/data/experiment/html/airway.html) package developped by **Michael Love**.  
It contains data from a RNA-Seq experiment on four human airway smooth muscle cell lines treated with dexamethasone. Details on the gene model and read counting procedure are provided in the package vignette.  
The experimental part and data generation was done by **Himes BE, Jiang X, Wagner P, Hu R, Wang Q, Klanderman B, Whitaker RM, Duan Q, Lasky-Su J, Nikolos C, Jester W, Johnson M, Panettieri R Jr, Tantisira KG, Weiss ST, Lu Q**: '*RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells.*' PLoS One. 2014 Jun 13;9(6):e99625. PMID: 24926665. GEO: GSE52778.


## Support
If you are having any trouble with the application, contact **victor.lebars@univ-rennes1.fr** or the dedicated [Github issue page](https://github.com/IGDRion/deseq2-viewer/issues).

## Authors and acknowledgment
Author: Victor Le Bars

## License
This project is licensed under the [MIT License](LICENSE) - see the [LICENSE](LICENSE) file for details.
