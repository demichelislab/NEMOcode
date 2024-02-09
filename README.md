
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Non-invasive detection of neuroendocrine prostate cancer through targeted cell-free DNA methylation

<!-- badges: start -->

[![License:
MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://cran.r-project.org/web/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8431083.svg)](https://doi.org/10.5281/zenodo.8431083)

<!-- badges: end -->

![NEMO](./img/schematic.png) *Schematic of NEMO panel design,
validation, and application*

## Description

This repository contains the required code to perform circulating Tumor
Content estimation (TC) and Phenotype Evidence score (PE score)
estimation from DNA methylation data generated (or masked) using the
NEMO assay. The NEMO assay is a tailored sequencing panel to monitor
CRPC through cell-free DNA (cfDNA) methylation and ultimately provide
detection of treatment-induced neuroendocrine prostate cancer. The
following steps estimate those clinically relevant quantities starting
from DNA methylation counts obtained from BS-seq techniques. For more
details, please see the original manuscript ([link](link)).

Notice that, in principle, masked whole genome bisulfite sequencing data
can be used to generate per-region Beta values. Even though the signal
quality will likely be lower due to limited coverage, we successfully
applied our method to a collection of WGBS and eRRBS datasets obtaining
coherent results.

Notably, the presented scripts are tailored for the NEMO panel design
and were not developed to be more general deconvolution or tumor
estimation strategies. While the same principles can be applied to other
tumor types or settings, our approach is not guaranteed to work beyond
its scope.

### Data availability

Preprocessed data, including DNA methylation, counts for samples
profiled with NEMO, and per-region DNA methylation values for all
samples analyzed, are available in Zenodo
[link](http://dx.doi.org/10.5281/zenodo.8431083). Raw WGBS data are
available from dbgap (WGBS, phs001752.v2-ctDNA_NE_CRPC_v2). Raw
sequencing data for samples profiled with NEMO are available upon
request due to privacy concerns, but anonymized bam files and processed
DNA methylation calls are available through GEO (GSE245745). All results
from the paper can be produced using only DNA methylation counts data.

A minimal guide to obtain NEMO estimation from raw data (DNA methylation
counts from Bismark) is outlined below.

### Data preparation

Raw data has been processed using a standard pipeline for DNA
methylation (nf-core/methylseq, version 1.6 (
[link](https://nf-co.re/methylseq/1.6/docs/usage) ). DNA methylation
counts in bismark format, and hg19 coordinates are required as a
starting input. Examples of such input files can be found in
`inst/extdata/`. Other types of inputs can be adapted, provided that per
CpG DNA methylation counts are available. The first functions
`collect_beta_panel()` and `collapse_stats()` use the input counts to
estimate a region-wise DNA methylation value (Beta). For each sample and
each region in the NEMO design, the average Beta value is computed as
the fraction of methylated CpGs over the total (first per CpG, then per
region). Our analyses retain only CpGs with at least 10x coverage, even
though this constraint can be relaxed. The output of this step can be
formatted into a matrix that recapitulates the average DNA methylation
level for each sample and each informative region.

![NEMO](./img/est_schematic.png) *Tumor content and Phenotype Evidence
score inference strategies*

### Inference of Tumor Content

53 regions of the NEMO panel have been selected to estimate tumor
content in circulation. Those regions are differentially methylated
between CRPC samples (no matter their class) and white blood cells. As
white blood cells are the major contributors to the healthy cfDNA pool,
those regions provide the ideal signal to estimate the relative
contribution of tumor derived DNA in circulation. The estimation
procedure leverages reference values from a collection of healthy cfDNA
samples (Fox-Fisher et al. eLIFE 2021) and the fact that each region has
opposite states in CRPC and healthy cfDNA. The function
`compute_score_ic()` produces a tumor content estimation for each sample
and a stability interval which results from the iterative downsampling
of the informative regions. Our tests indicate that this strategy can
quantify tumor presence up to 3%.

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
healthy cfDNA signal, as it is equal to one minus the tumor content,
while a non-informative prior is used for the NE and Adeno components.
The relative error of the estimation is computed as the maximum standard
error across NE and Adeno fractions, divided by the total tumor content
(currently, this statististic is not used). Notice that the PE score is
not reliable when the tumor content is below 3%, while we recommend
caution in interpreting the PE score between 3% and 10% of tumor
content.

### Installation

You can install the current version of NEMOcode from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
install.packages("brms")
devtools::install_github("demichelislab/NEMOcode")
```

### Requirements

These package has been tested with R version `4.1.2` and Ubuntu 22.04. A
series of additioal packages will be installed to make `NEMOcode` work.
Those dependencies are listed in the `DESCRIPTION` file. The `brms`
package might require additional care for the installation, see this
[link](https://learnb4ss.github.io/learnB4SS/articles/install-brms.html)
for help.  
Installation of the package should take few minutes.

Notice that MacOS systems might require the additional installation of
`Xcode`, `gfortran`, `StanHeaders`, and R dev binaries using
`macrtools`.

## Example

``` r
library(dplyr)
library(tidyr)
library(NEMOcode)
rstan::rstan_options(javascript=FALSE)
```

### Preparation

First, we load a small set of data included in the package. Those
samples have been profiled with the NEMO assay, and this data is in a
standard form of raw counts. The three examples include a CRPC organoid
with stem-cell like features (`PM155_P`), a cell-free DNA sample from an
healthy donor (`HD1`), a neuroendocrine cell line (`NCI-H660`), and the
VCaP CRPC cell line (`VCaP`).

``` r

fpath1 <- system.file("extdata", "PM155_P-run_pilotCells_unif.RData", package="NEMOcode")
fpath2 <- system.file("extdata", "HD1-run_pilot_cfDNA_unif.RData", package="NEMOcode")
fpath3 <- system.file("extdata", "VCaP-run_pilotCells_unif.RData", package="NEMOcode")
fpath4 <- system.file("extdata", "NCI-H660-run_pilotCells2_unif.RData", package="NEMOcode")

exfile1 = local(get(load(fpath1)))
exfile2 = local(get(load(fpath2)))
exfile3 = local(get(load(fpath3)))
exfile4 = local(get(load(fpath4)))

head(exfile1)
#>    chr   start     end      beta meth unmeth
#> 1 chr1 2036143 2036143 0.4210526    8     11
#> 2 chr1 2036162 2036162 0.5454545   12     10
#> 3 chr1 2036166 2036166 0.3636364    8     14
#> 4 chr1 2036185 2036185 0.4400000   11     14
#> 5 chr1 2036190 2036190 0.4166667   10     14
#> 6 chr1 2036202 2036202 0.2500000    6     18
```

The first step of the analysis is to reduce the input data to region
level DNA methylation, using the panel design as reference. For each
region, the average DNA methylation (Beta) is computed. The high local
autocorrelation of DNA methylation supports the idea of collapsing
contiguous CpGs to obtain a robust signal.

``` r

list_counts = list("PM155_P" = exfile1, "HD1" = exfile2, "VCaP" = exfile3, "NCI-H660" = exfile4)
pdes = NEMOcode::panel_design
coll_regs = list()

for (i in names(list_counts)){
    
    colldes = collect_beta_panel(cov_file = list_counts[[i]], panel_design = pdes)
    statsdes = collapse_stats(colldes, ssid = i)
    coll_regs[[i]] = statsdes
}

coll_regs = bind_rows(coll_regs) %>% arrange(sampleID, reg_id)
head(coll_regs)
#> # A tibble: 6 × 8
#> # Groups:   region_set [2]
#>   region_set           reg_id mean_meth sd_meth n_cpgs tot_cov tot_meth sampleID
#>   <chr>                <chr>      <dbl>   <dbl>  <int>   <dbl>    <dbl> <chr>   
#> 1 rockermeth_NEvsAdeno chr10…    0.210  NA           1     186       39 HD1     
#> 2 NE_Adeno             chr10…    0.0520  0.0251     21    4697      242 HD1     
#> 3 rockermeth_NEvsAdeno chr10…    0.96   NA           1     175      168 HD1     
#> 4 rockermeth_NEvsAdeno chr10…    0.974  NA           1     154      150 HD1     
#> 5 rockermeth_NEvsAdeno chr10…    0.981  NA           1     157      154 HD1     
#> 6 rockermeth_NEvsAdeno chr10…    0.975  NA           1     200      195 HD1
```

Now we can turn the result into matrix form, using the average Beta
value per region as entry. This will be the main input for tumor content
estimation and PE score estimation.

``` r

test_mat = coll_regs %>% 
    ungroup() %>% 
    select(sampleID, mean_meth, reg_id) %>% 
    pivot_wider(names_from = sampleID, values_from = mean_meth) %>% 
    as.data.frame()

to_matrix<-function(x) {
    n = x[,1]
    m<-as.matrix(x[,2:ncol(x)])
    rownames(m)<-x[,1]
    m
}


test_mat = to_matrix(test_mat)
head(test_mat)
#>                                  HD1   NCI-H660   PM155_P      VCaP
#> chr10:119590052-119590055 0.20967742 0.93571429 0.9687500 0.8888889
#> chr10:123923781-123924185 0.05196424 0.04188003 0.9762074 0.9724014
#> chr10:134562610-134562613 0.96000000 0.96153846 0.9117647 0.9677419
#> chr10:134563052-134563055 0.97402597 0.96385542 1.0000000 0.9090909
#> chr10:134563060-134563064 0.98089172 0.96590909 0.9714286 1.0000000
#> chr10:134563098-134563102 0.97500000 0.96969697 0.9736842 0.9459459
```

### Tumor content estimation

The next step is tumor content estimation. In this case, only the tumor
informative regions are used (`mCRPC_PBMC`, as they have CRPC specific
signal against the white blood cell background).

``` r

ptc = pdes %>% filter(region_set == "mCRPC_PBMC")
control_mat = NEMOcode::control_mat

test_mat_tc = test_mat[ptc$reg_id, ]
control_mat = control_mat[ptc$reg_id, ]

state = ptc %>% select(reg_id, meth_state)

res = compute_ci_confidence(
    tumor_mat = test_mat_tc, 
    control_mat = control_mat,
    reg_df = state,
    nsub = 100,
    quant_prob = 0.05,
    frac_sub = 0.5
    )

res
#> # A tibble: 4 × 9
#>   SampleName est_mu  est_sd est_min est_max q025_tc q975_tc ci_lower ci_upper
#>   <chr>       <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>    <dbl>    <dbl>
#> 1 HD1             0 0         0           0       0       0        0        0
#> 2 NCI-H660        1 0.00125   0.991       1       1       1        1        1
#> 3 PM155_P         1 0         1           1       1       1        1        1
#> 4 VCaP            1 0         1           1       1       1        1        1
```

As expected, the cfDNA sample from an healthy donor is estimated to have
a tumor content (`est_mu`) below 3%, our sensitivity threshold, and thus
it can be considered negative for tumor. On the other hand, all CRPC
models produce a tumor content of 1, as expected.

### Phenotype Evidence score estimation

Once the tumor content is known, the final step is to estimate the
relative contributions of CRPC-NE, CRPC-Adeno, and cfDNA. In this case,
another set of informative regions is used (a collection of regions
resulted from differential methylation analysis between CRPC-Adeno and
CRPC-NE samples). Importantly, the tumor content is known and it is fed
to the function, informing the expected relative fraction of CRPC signal
over the total. In the `compute_all()` function it is possible to
indicate a set of samples for which `cfDNA` is the expected background,
instead of `PBMC` (default). Furthermore, one can enforce the tumor
content estimation to be 1 by adding the sample name to `atlas_samples`.

``` r

ppe = pdes %>%
    filter(
        region_set %in% c(
            "rockermeth_NEvsAdeno",
            "NE_Adeno",
            "CpG_site_DNAmeth_classifierJCI_2020"
        )
    )

refDist = NEMOcode::refDist
rownames(refDist) = refDist$reg_id

## Use only common rows and ensure that the order is the same
kk = base::intersect(rownames(refDist), ppe$reg_id)
refDist = refDist[kk,]
obs_all = test_mat[kk, ]

result = compute_all(
    obs_all,
    ref_mat = refDist,
    atlas_tc = res,
    atlas_samples = c(),
    cfdna_ids = c("HD1"),
    sequential = T
)

result
#>         immune     adeno        ne rel_error var_score SampleName       pes
#> 1 1.0000000000        NA        NA        NA        NA        HD1        NA
#> 2 0.0003966096 0.2369646 0.7626388 0.3647717 0.1508964   NCI-H660 0.7629414
#> 3 0.0023765074 0.4889418 0.5086817 0.3617121 0.1398005    PM155_P 0.5098934
#> 4 0.0000000000 0.8854295 0.1145705 0.2533096 0.1038627       VCaP 0.1145705
#>      pes_lw    pes_up tc_est tc_lw tc_up quality_flag
#> 1        NA        NA      0     0     0        FALSE
#> 2 0.7629414 0.7629414      1     1     1         TRUE
#> 3 0.5098934 0.5098934      1     1     1         TRUE
#> 4 0.1145705 0.1145705      1     1     1         TRUE
```

The results report the relative contributions of the three expected
populations `ne`, `adeno` and `immune`. The `pe` score is computed as
`ne/(ne+adeno)`, and reflects the degree of neuroendocrine
differentiation of tumor derived DNA.

Notice that the runtime for each sample can be in the order of few
minutes, as the bayesian regression employs numerical methods to
estimate regression parameters. For this reason, the analysis of a large
dataset might require parallelization over multiple cores on a machine
with good computational power (parallelization is already implemented in
`compute_all()` function. `sequential=TRUE` sets sequential computation.
`nclust` can be used to indicate the number of threads to use).

The relative error and variability score are simple metrics to measure
reconstruction accuracy. Currently, they are not used but were
considered to set the lower limit of PE score estimation. Their usage
will possibly be considered in future implementation. The same holds for
the confidence interval of the PE score, which is currently a
placeholder. The quality flag variable simply reports weather the sample
has the minimum tumor content to perform an estimation (3%) or not.

To reproduce the score reported in the paper simply run the scoring
functions on the corresponding data deposited in
[Zenodo](https://doi.org/10.5281/zenodo.7887998). Computed scores for
analyzed samples are also available as Supplementary Data in our study
(Franceschini et al., …).

## Docker

A singularity image containing `NEMOcode` and the required packages can
be build using the following recipe:

``` singularity
Bootstrap: docker
From: ubuntu:22.04
IncludeCmd: yes

%environment
R_VERSION=4.1.2
export R_VERSION
R_CONFIG_DIR=/etc/R/
export R_CONFIG_DIR
export LC_ALL=C
export PATH=$PATH

%labels
Author Francesco Orlando
Version v1.0.0
R_Version 4.2.2
build_date 08/02/2024
R_bioconductor True

%apprun R
exec R "$@"

%apprun Rscript
exec Rscript "$@"

%runscript
exec R "$@"

%post
  apt-get update
  apt-get install -y apt-transport-https apt-utils software-properties-common
  apt-get install -y add-apt-key
  export DEBIAN_FRONTEND=noninteractive
  ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime
  apt-get install -y tzdata
  dpkg-reconfigure --frontend noninteractive tzdata

  apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
  add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
  apt-get update

  apt-get install -y gcc fort77 aptitude
  aptitude install -y g++
  aptitude install -y gfortran
  apt-get install -y r-base r-base-dev

  apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    make \
    g++ \
    libboost-all-dev \
    libglpk-dev \
    libgmp-dev \
    libxml2-dev
  echo 'CXX14 = g++ -std=c++1y -Wno-unused-variable -Wno-unused-function -fPIC' >> /etc/R/Makeconf
  
  # Install the tidyverse, devtools, rstan, and brms packages in R
R -e "options(repos = list(CRAN = 'https://cloud.r-project.org/')); \
         install.packages(c('tidyverse', 'rstan', 'BiocManager', 'brms', 'remotes'), dependencies=TRUE)"
R -e "remotes::install_github('DesiQuintans/librarian')"
R -e "install.packages('igraph', dependencies = TRUE)"
R -e "BiocManager::install(update=FALSE)"
R -e "BiocManager::install('Biobase', update=FALSE)"
R -e "librarian::shelf(GenomicRanges, rtracklayer, S4Vectors)"
R -e "remotes::install_github('demichelislab/NEMOcode', dependencies = TRUE)"
```
