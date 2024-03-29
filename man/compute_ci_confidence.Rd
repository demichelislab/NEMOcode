% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deconv.R
\name{compute_ci_confidence}
\alias{compute_ci_confidence}
\title{Compute the fraction of unexpected epithelial cell-free DNA in a sample with subsampling confidence interval}
\usage{
compute_ci_confidence(
  tumor_mat,
  control_mat,
  reg_df,
  nsub = 100,
  quant_prob = 0.05,
  frac_sub = 0.5
)
}
\arguments{
\item{tumor_mat}{Matrix of plasma samples in which estimate the epithelial (tumoral) fraction, regions on rows, samples on columns}

\item{control_mat}{Matrix of PBMC/cfDNA profiles, regions on rows, samples on columns}

\item{reg_df}{Dataframe of informative regions with 2 columns reg_id and meth state (-1/+1) for hypo and hyper respectively, reporting the differential between CRPC samples and PBMC/cfDNA}

\item{nsub}{number of subsampling iterations to estimate stability (default = 100)}

\item{quant_prob}{the quantile of probabilities to produce stability interval via subsampling (0.05 = 95\% CI)}

\item{frac_sub}{the fraction of informative regions to use for subsampling iterations}
}
\value{
A dataframe with tumor content estimation
}
\description{
Compute the fraction of unexpected epithelial cell-free DNA in a sample with subsampling confidence interval
}
