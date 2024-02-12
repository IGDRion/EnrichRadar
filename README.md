# DESeq2 viewer

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Description
This R shiny application allows to view the output of DESeq2 and make interactive filters over the data in an easy way. It is also able to produce an enrichment analysis via Gprofiler2 and return an interactive figure.

## Installation

- Git clone the repository:

```bash
git clone https://github.com/IGDRion/deseq2-viewer.git
```

- Install the required R libraries using e.g. with `conda`

```bash
conda install conda-forge::r-shiny r::r-dt conda-forge::r-ggplot2  plotly::plotly conda-forge::r-gprofiler2 r-dipsaus conda-forge::r-dplyr conda-forge::r-shinyjs
```

- Launch the app:

```bash
Rscript /path/to/app.R
```

It should open the application on your web browser : `http://127.0.0.1:7449`.

## Usage
To load a DESeq2 output file, click on the browse button and select the file you would like to view.
Once the file is loaded, you can filter your data :
- Choose wich threshold you want to use for Log2FoldChange and p-adujsted to filter your data
- Select which type of differential expression you want to see : "Both" shows down and up regulated genes, "Down" shows only down regulated genes, and "Up" shows only up regulated genes
- Two checkboxes allows you to filter if you want to see only protein coding genes & only known genes
- "Save current selection" button allows to save the dataframe with all the modification done to it (sorting of rows are not saved however)
- Gprofiler button generate an interactive graph with gprofiler2 (tool for multi-databases enrichment analysis) 

A Toy example is available to try the application : **ToyExample.csv**
The data are from the [**airway** dataset](https://bioconductor.org/packages/release/data/experiment/html/airway.html) package developped by **Michael Love**.  
This data are from a RNA-Seq experiment on four human airway smooth muscle cell lines treated with dexamethasone. Details on the gene model and read counting procedure are provided in the package vignette.  
The experiment was done by **Himes BE, Jiang X, Wagner P, Hu R, Wang Q, Klanderman B, Whitaker RM, Duan Q, Lasky-Su J, Nikolos C, Jester W, Johnson M, Panettieri R Jr, Tantisira KG, Weiss ST, Lu Q**: '*RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells.*' PLoS One. 2014 Jun 13;9(6):e99625. PMID: 24926665. GEO: GSE52778.


## Support
If you are having any trouble with the application, contact **victor.lebars@univ-rennes1.fr** or the dedicated [Github issue page](https://github.com/IGDRion/deseq2-viewer/issues).

## Authors and acknowledgment
Author: Victor Le Bars

## License
This project is licensed under the [MIT License](LICENSE) - see the [LICENSE](LICENSE) file for details.
