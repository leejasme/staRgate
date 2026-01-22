# staRgate
A density-based automated gating pipeline for flow cytometry data to characterize the lineage, differentiation, and functional states of T-cells

This is the GitHub repo for the _staRgate_ R package. 

# Authors

Jasme Lee, Matthew Adamow, Colleen Maher, Xiyu Peng, Phillip Wong, Fiona Ehrich, Michael A. Postow, Margaret K. Callahan, Ronglai Shen, Katherine S. Panageas

# Installation



**Bioconductor version**

+ [Please refer to Bioconductor for full details on installation guidelines](https://www.bioconductor.org/install/)

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("staRgate")
```

**GitHub version (developmental)**

```
devtools::install_github("leejasme/staRgate")
```

# Tutorial 

A full example on [how to run the pipeline is on the webpage](https://leejasme.github.io/staRgate/articles/vignette_run_pipeline.html)


# Funding Sources

This research was funded by the V foundation, MSK-MIND, NIH R01CA276286, and NIH P30CA008748.

# Contact

Jasme Lee (leej22@mskcc.org)

