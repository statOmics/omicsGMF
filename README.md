# omicsGMF

omicsGMF is a Bioconductor package that uses the sgdGMF-framework of the sgdGMF package for highly performant 
    and fast matrix factorization that can be used for dimensionality reduction, visualization and imputation of omics
    data. It considers data from the general exponential family as input, and therefore suits the use of both RNA-seq
    (Poisson or Negative Binomial data) and proteomics data (Gaussian data). omicsGMF allows for:

- Dimensionality reduction directly on (single-cell) RNA-sequencing count data considering a Poisson or negative binomial family. Therefore, no prior normalization is needed as with conventional PCA.
- Dimensionality reduction when dealing with missing data. This makes it suitable for (single-cell) proteomics data when considering a Gaussian family.
- Correction of sample-level (e.g. batch effects) and feature-level covariates (e.g. GC-content). Therefore, the dimensionality reduction obtained is useful when visualizing or clustering cells in which technical variability is present.
- Imputation of missing values by the estimated means, which also take into account the known covariates.

For more information, be sure to check out:
- bioRxiv preprint on omicsGMF applied to proteomics data for dimensionality reduction and imputation of missing values:
- arXiv preprint on the technical manuscript of sgdGMF framework, applied to single-cell RNA-sequencing: https://arxiv.org/abs/2412.20509



## Installation instructions

To install the development version, run;

```{r 'install_dev', eval = FALSE}
devtools::install_github("statOmics/omicsGMF")
```

The installation should only take a few seconds.
The dependencies of the package are listed in the DESCRIPTION file of the package.

omicsGMF is submitted to Bioconductor and should be soon available. 

## Issues and bug reports

Please use https://github.com/statOmics/omicsGMF/issues to submit issues, bug reports, and comments.


