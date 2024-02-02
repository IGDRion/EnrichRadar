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
conda install conda-forge::r-shiny r::r-dt conda-forge::r-ggplot2  plotly::plotly conda-forge::r-gprofiler2 r-dipsaus conda-forge::r-dplyr
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
The data of for this example were extracted from [this analysis](https://www.ebi.ac.uk/gxa/experiments/E-GEOD-50760/Supplementary%20Information). It is a RNA-seq  analysis of 54 samples (normal colon, primary tumor, and liver metastases) from 18 colorectal cancer patients done by *Kim SK, Kim SY, Kim JH, Roh SA, Cho DH et al. (2014)* [A nineteen gene-based risk score classifier predicts prognosis of colorectal cancer patients](https://europepmc.org/abstract/MED/25049118). Out of this dataset, only colon and colon epithelium samples were kept, and their expression profile were compared using DESeq2 to produce the **ToyExample** csv dataset.  
Please note that this dataset is not intended to be representative of a coherent analysis, but rather to be used to test the R shiny application.


## Support
If you are having any trouble with the application, contact **victor.lebars@univ-rennes1.fr** or the dedicated [Github issue page](https://github.com/IGDRion/deseq2-viewer/issues).

## Authors and acknowledgment
Author: Victor Le Bars

## License
This project is licensed under the [MIT License](LICENSE) - see the [LICENSE](LICENSE) file for details.
