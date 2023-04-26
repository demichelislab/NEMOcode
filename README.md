
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Companion code for the NEMO sequencing panel

<!-- badges: start -->

[![codecov](https://codecov.io/gh/GMFranceschini/NEMOcode/branch/master/graph/badge.svg?token=M6G0L06263)](https://codecov.io/gh/GMFranceschini/NEMOcode)
<!-- badges: end -->

## Installation

You can install the current version of NEMOcode from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("GMFranceschini/NEMOcode")
```

## Description

### NEMO tumor content and PE score estimation

This repository contains the required code to compute Tumor Content
estimation (TC) and Phenotype Evidence score (PE score) estimation from
DNA methylation data generated or masked using the NEMO-panel design.
The NEMO-panel is a tailored sequencing assay to monitor CRPC through
cell-free DNA (cfDNA) methylation and ultimately provide detection of
treatment-induced neuroendocrine prostate cancer. Starting from DNA
methylation counts obtained from BS-seq techniques, the following steps
estimate those clinically relevant quantities. For more details, please
see the original manuscript [link](link).

Notice that, in principle, masked whole genome bisulfite sequencing data
can be used to generate per-region Beta values, even though the signal
quality will likely be lower due to limited coverage.

Importantly, the presented scripts are tailored for the NEMO panel
design and where not developed to be more general deconvolution or tumor
estimation strategies.

### Requirements

These scripts have been tested with R version `4.1.2`.

### Data availability

Preprocessed data are available in Zenodo [link](link). Raw data are
available from GEO (NEMO panel) and SRA (WGBS,
phs001752.v2-ctDNA_NE_CRPC_v2).

A minimal guide to obtain NEMO estimation from raw data is outlined
below.

### Data preparation

Data has been processed using a standard pipeline for DNA methylation
(`nf-core/methylseq`, version 1.6 [link](link)). DNA methylation counts
in bismark format and hg19 coordinates are required as a starting input.
Other types of inputs can be adapted, provided that per CpG DNA
methylation counts are available. Examples of such input files can be
found in `inst/dataext/rawcov`. The first function `collapse_beta()`
uses the input counts to estimate a region-wise DNA methylation value
(Beta). For each sample and each region in the NEMO design, the average
Beta value is computed as the fraction of methylated CpGs over the total
(first per CpG, then per region), keeping only CpG sites with coverage
greater than 10. The result of this step is a matrix with samples in the
columns and NEMO regions in rows.

### Inference of Tumor Content

53 regions of the NEMO panel have been selected to estimate tumor
content in circulation. Those regions are differentially methylated
between CRPC samples (independently from their phenotype) and white
blood cells. As white blood cells are the major contributors to the
healthy cfDNA pool, those regions provide the ideal signal to estimate
the relative contribution of tumor derived DNA in circulation. The
estimation procedure leverages reference values from a collection of
healthy cfDNA samples (Fox-Fisher et al. eLIFE 2021) and the fact that
each region has opposite states in CRPC and healthy cfDNA. The function
`compute_score_ic()` produces a tumor content estimation for each sample
and a stability interval which results from the iterative downsampling
of the informative regions.

### Inference of PE score

Once the tumor content is known, the main goal of the analysis is to
infer the presence of CRPC-NE derived DNA. As differential DNA
methylation between CRPC-Adeno and CRPC-NE has been previously
characterized, we leverage a set of informative regions to achieve this
goal. The inference procedure builds upon the tumor content estimation
and aims at reconstructing the observed signal as a linear combination
of CRPC-NE, CRPC-Adeno, and healthy cfDNA signal. To this end, an atlas
with reference profiles generated with high tumor-content tissue samples
has been generated and is available in the `resource` folder. The
reconstruction task is performed through a Bayesian linear regression
`compute_all()`/`compute_evidence_brms()`, which computes the most
likely fractions of each cell type. A strong prior is used for the
healthy cfDNA signal, as it is equal to $(1 - TC)$, while a
non-informative prior is used for the NE and Adeno components. The
relative error of the estimation is computed as the maximum standerd
error across NE and Adeno fractions, divided by the total tumor content.
The PE score is not reliable when the tumor content is below 5%, while
we recommend caution in interpreting the PE score between 5% and 10% of
TC.

``` r
library(NEMOcode)
## basic example code
```
