% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/omicsGMF-package.R
\docType{package}
\name{omicsGMF-package}
\alias{omicsGMF}
\alias{omicsGMF-package}
\title{omicsGMF: Dimensionality reduction of (single-cell) omics data in R using omicsGMF}
\description{
omicsGMF is a Bioconductor package that uses the sgdGMF-framework of the \code{sgdGMF} package for highly performant and fast matrix factorization that can be used for dimensionality reduction, visualization and imputation of omics data. It considers data from the general exponential family as input, and therefore suits the use of both RNA-seq (Poisson or Negative Binomial data) and proteomics data (Gaussian data). It does not require prior transformation of counts to the log-scale, because it rather optimizes the deviances from the data family specified. Also, it allows to correct for known sample-level and feature-level covariates, therefore enabling visualization and dimensionality reduction upon batch correction. Last but not least, it deals with missing values, and allows to impute these after matrix factorization, useful for proteomics data. This Bioconductor package allows input of SummarizedExperiment, SingleCellExperiment, and QFeature classes.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/statOmics/omicsGMF}
  \item Report bugs at \url{https://github.com/statOmics/omicsGMF/issues}
}

}
\author{
\strong{Maintainer}: Alexandre Segers \email{alexandresegers@outlook.com}

Other contributors:
\itemize{
  \item Cristian Castiglione \email{cristian.castiglione@unipd.it} [contributor]
  \item Christophe Vanderaa \email{christophe.vanderaa@ugent.be} [contributor]
  \item Davide Risso \email{davide.risso@unipd.it} [contributor]
  \item Lieven Clement \email{lieven.clement@ugent.be} [contributor]
}

}
\keyword{internal}
