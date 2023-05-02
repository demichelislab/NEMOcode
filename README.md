
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Non-invasive detection of neuroendocrine prostate cancer through targeted cell-free DNA methylation

<!-- badges: start -->

[![codecov](https://codecov.io/gh/GMFranceschini/NEMOcode/branch/master/graph/badge.svg?token=M6G0L06263)](https://codecov.io/gh/GMFranceschini/NEMOcode)
[![License:
MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://cran.r-project.org/web/licenses/MIT)
[![doi](https://img.shields.io/badge/doi-000000-green.svg)](https://)
<!-- badges: end -->

![NEMO](./img/schematic.png)

## Description

### NEMO tumor content and PE score estimation

This repository contains the required code to compute circulating Tumor
Content estimation (TC) and Phenotype Evidence score (PE score)
estimation from DNA methylation data generated or masked using the
NEMO-panel design. The NEMO-panel is a tailored sequencing assay to
monitor CRPC through cell-free DNA (cfDNA) methylation and ultimately
provide detection of treatment-induced neuroendocrine prostate cancer.
Starting from DNA methylation counts obtained from BS-seq techniques,
the following steps estimate those clinically relevant quantities. For
more details, please see the original manuscript [link](link).

Notice that, in principle, masked whole genome bisulfite sequencing data
can be used to generate per-region Beta values. Even though the signal
quality will likely be lower due to limited coverage, we successfully
applied our method to a collection of WGBS and eRRBS datasets.

Importantly, the presented scripts are tailored for the NEMO panel
design and where not developed to be more general deconvolution or tumor
estimation strategies. While the same principles can be applied to other
tumor types, there is no guarantee that our approach will work beyond
its scope.

### Data availability

Preprocessed data are available in Zenodo [link](link). Raw data are
available from GEO (NEMO panel) and SRA (WGBS,
phs001752.v2-ctDNA_NE_CRPC_v2).

A minimal guide to obtain NEMO estimation from raw data (DNA methylation
counts from Bismark) is outlined below.

### Data preparation

Raw data has been processed using a standard pipeline for DNA
methylation (`nf-core/methylseq`, version 1.6 [link](link)). DNA
methylation counts in bismark format and hg19 coordinates are required
as a starting input. Examples of such input files can be found in
`inst/extdata/`. Other types of inputs can be adapted, provided that per
CpG DNA methylation counts are available. The first functions
`collect_beta_panel` and `collapse_stats` use the input counts to
estimate a region-wise DNA methylation value (Beta). For each sample and
each region in the NEMO design, the average Beta value is computed as
the fraction of methylated CpGs over the total (first per CpG, then per
region). In our analyses, we retain only CpGs with at least 10x
coverage, even though this constraint can be relaxed. The output of this
step can be formatted into a matrix that recapitulates the average DNA
methylation level for each sample and each informative region.

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
healthy cfDNA signal, as it is equal to $(1 - TC)$, while a
non-informative prior is used for the NE and Adeno components. The
relative error of the estimation is computed as the maximum standerd
error across NE and Adeno fractions, divided by the total tumor content.
The PE score is not reliable when the tumor content is below 3%, while
we recommend caution in interpreting the PE score between 3% and 10% of
TC.

### Installation

You can install the current version of NEMOcode from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("GMFranceschini/NEMOcode")
```

### Requirements

These scripts have been tested with R version `4.1.2` and Ubuntu 22.04.
A series of packages will be installed to make `NEMOcode` work. Those
dependencies are listed in the `DESCRIPTION` file. The `brms` package
might require additional care for the installation, see this
[link](https://learnb4ss.github.io/learnB4SS/articles/install-brms.html)
for help. Installation of the package should take few minutes.

``` r
library(tidyverse)
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.1.2     ✔ readr     2.1.4
#> ✔ forcats   1.0.0     ✔ stringr   1.5.0
#> ✔ ggplot2   3.4.2     ✔ tibble    3.2.1
#> ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
#> ✔ purrr     1.0.1     
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
library(NEMOcode)
## basic example code
```

First, we load a small set of data included in the package. Those
samples have been profiled with the NEMO assay, and this data is in a
standard form of raw counts.

``` r

fpath1 <- system.file("extdata", "PM155_P-run_pilotCells_unif.RData", package="NEMOcode")
fpath2 <- system.file("extdata", "HD1-run_pilot_cfDNA_unif.RData", package="NEMOcode")
fpath3 <- system.file("extdata", "VCaP-run_pilotCells_unif.RData", package="NEMOcode")

exfile1 = local(get(load(fpath1)))
exfile2 = local(get(load(fpath2)))
exfile3 = local(get(load(fpath3)))

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
autocorrelation supports the idea of collapsing contiguous CpGs to
obtain a robust signal.

``` r

list_counts = list("PM155_P" = exfile1, "HD1" = exfile2, "VCaP" = exfile3)
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
value per region as entry.

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
#>                                  HD1   PM155_P      VCaP
#> chr10:119590052-119590055 0.20967742 0.9687500 0.8888889
#> chr10:123923781-123924185 0.05196424 0.9762074 0.9724014
#> chr10:134562610-134562613 0.96000000 0.9117647 0.9677419
#> chr10:134563052-134563055 0.97402597 1.0000000 0.9090909
#> chr10:134563060-134563064 0.98089172 0.9714286 1.0000000
#> chr10:134563098-134563102 0.97500000 0.9736842 0.9459459
```

The next step is tumor content estimation. In this case, only the TC
informative regions are used (`mCRPC_PBMC"`, as they have CRPC specific
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
#> # A tibble: 3 × 9
#>   SampleName est_mu  est_sd est_min est_max q025_tc q975_tc ci_lower ci_upper
#>   <chr>       <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>    <dbl>    <dbl>
#> 1 HD1        0.0279 0.00671  0.0200  0.0374  0.0209  0.0365   0.0209   0.0365
#> 2 PM155_P    1      0        1       1       1       1        1        1     
#> 3 VCaP       1      0        1       1       1       1        1        1
```

Once the tumor content is known, the final step is to estimate the
relative contributions of CRPC-NE, CRPC-Adeno, and cfDNA. In this case,
another set of informative regions is used (a collection of regions
resulted from differential methylation analysis between CRPC-Adeno and
CRPC-NE samples). Importantly, the tumor content is known and it is fed
to the function, informing the expected relative fraction of CRPC signal
over the total.

``` r

ppe = pdes %>% 
    filter(region_set %in% c("rockermeth_NEvsAdeno", "NE_Adeno", "CpG_site_DNAmeth_classifierJCI_2020"))

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
    cfdna_ids = c(),
    sequential = T
)

result
#>        immune       adeno         ne  rel_error  var_score SampleName       pes
#> 1 0.973586277 0.001861452 0.02455227 0.08750612 0.04714514        HD1        NA
#> 2 0.001911564 0.487815726 0.51027271 0.36210107 0.12367707    PM155_P 0.5112500
#> 3 0.000000000 0.885031668 0.11496833 0.25340946 0.10335912       VCaP 0.1149683
#>      pes_lw    pes_up     tc_est      tc_lw      tc_up quality_flag
#> 1        NA        NA 0.02790538 0.02086163 0.03647445        FALSE
#> 2 0.5112500 0.5112500 1.00000000 1.00000000 1.00000000         TRUE
#> 3 0.1149683 0.1149683 1.00000000 1.00000000 1.00000000         TRUE
```

Notice that the runtime for each sample can be in the order of few
minutes, as the bayesian regression employs numerical methods to
estimate regression parameters. For this reason, the analysis of a large
dataset might require parallelization over multiple cores on a machine
with good computational power (parallelization is already implemented in
`compute_all()` function).

To reproduce the score reported in the paper simply run the scoring
functions on the corresponding data deposited in Zenodo
[Zenodo](https://github.com/).
