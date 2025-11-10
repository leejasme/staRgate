# staRgate

A density-based automated gating pipeline for flow cytometry data to
characterize the lineage, differentiation, and functional states of
T-cells

This GitHub stores the {staRgate} R package.

# Installation

The {staRgate} package relies on a few Biocondcutor R packages. Before
installing {staRgate}, first setup Bioconductor and install all
packages.

- {flowCore} and {flowWorkspace} are dependencies for the {staRgate}
  package
- {openCyto}, {flowAI}, {ggcyto} are not required to run the functions
  of {staRgate} but are used in the full gating pipeline as shown in the
  [Tutorial](https://leejasme.github.io/staRgate/articles/vignette_run_pipeline.html)

[Please refer to Bioconductor for full details on installation
guidelines](https://www.bioconductor.org/install/)

    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

To install staRgate (currently install from GitHub):

    devtools::install_github("leejasme/staRgate")

# Tutorial

A full example on [how to run the pipeline is on the
webpage](https://leejasme.github.io/staRgate/articles/vignette_run_pipeline.html)

# Contact

Jasme Lee (<leej22@mskcc.org>)
