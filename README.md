<!-- README.md is generated from README.Rmd. Please edit that file -->

scEC
====

<!-- badges: start -->
<!-- badges: end -->

Single-cell Entropic Clustering (scEC) provides an entropy based
approach to normalisation, feature selection, differential expression
analysis and unsupervised clustering of single-cell RNA-sequencing data.

Installation
------------

You can install the the development version from
[GitHub](https://github.com/) with:

    # install.packages("devtools")
    devtools::install_github("mjcasy/scEC")

Workflow
--------

The basic workflow is demonstrated on the Tian et al (2018) single-cell
mixology data set, a mixture of three cancerous cell lines
(<a href="https://github.com/LuyiTian/sc_mixology" class="uri">https://github.com/LuyiTian/sc_mixology</a>).

    library(scEC)
    library(Matrix)

    load(paste0(Path, "CountsMatrix"))

    CountsMatrix <- CountsMatrix[rowSums(CountsMatrix) >= 100,]

    PopHet <- Population(CountsMatrix)
    Mean <- log10(rowMeans(CountsMatrix))
    N <- ncol(CountsMatrix)

    plot(Mean, PopHet, ylim = c(0, log(N)))

<img src="man/figures/README-example-1.png" width="100%" />

    GOI <- FeatureSelection(CountsMatrix)
    CountsMatrix <- CountsMatrix[GOI,]
    Mean <- Mean[GOI]

    Clustering <- Cluster(CountsMatrix, numClus = 3)
    InterType <- DifferentialExpression(CountsMatrix)

    plot(Mean, InterType)

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />
