#' Impute missing values based on the results of runSGD.
#'
#' @param x a numeric matrix of expression counts or mass spectrometry intensities containing missing values and with
#' features in the rows and samples in columns.
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such a matrix.
#' @param sgdGMF_reducedDims the output obtained by \code{'runSGD'} or \code{'calculateSGD'}. If \code{x} is a
#' \linkS4class{SingleCellExperiment}, \code{sgdGMF_reducedDims} is taken from \code{\link{reducedDim}(x, 'SGD')}.
#' @param assay.type Integer scalar or string indicating which assay of \code{x} contains the values of interest.
#' @param exprs_values Alias to \code{assay.type}.
#'
#' @details
#' Imputation is only possible after running \code{'runSGD'} using all features.
#' Therefore, \code{subset_row} or \code{ntop} should be set to NULL.
#'
#' @return
#' A SingleCellExperiment object is returned with an extra assay containing the imputed values.
#'
#' @name SGDImpute
#' @seealso
#' \code{\link[BiocSGDGMF]{runSGD}}, to conveniently obtain the matrix factorization.
#' \code{\link[BiocSGDGMF]{plotSGD}}, to conveniently visualize the results.
#'
#' @author Alexandre Segers
#'
#' @examples
#' example_sce <- mockSCE()
#'
#'
#' example_sce <- logNormCounts(example_sce)
#'
#' example_sce <- runSGD(example_sce)
#' reducedDimNames(example_sce)
#' head(reducedDim(example_sce))
NULL

.imputeMissingValues <- function(x, sgdGMF_reducedDims)
{
  means <- attr(sgdGMF_reducedDims,'B') %*% t(attr(sgdGMF_reducedDims,'X')) +
    attr(sgdGMF_reducedDims,'Z') %*% t(attr(sgdGMF_reducedDims,'Gamma')) +
    attr(sgdGMF_reducedDims, 'rotation') %*% t(sgdGMF_reducedDims)

  x[is.na(x)] <- means[is.na(x)]
  x
}


#' @export
#' @rdname SGDImpute
setMethod("SGDImpute", "ANY", .imputeMissingValues)

#' @export
#' @rdname SGDImpute
#' @importFrom SummarizedExperiment assay
setMethod("SGDImpute", "SummarizedExperiment", function(x, sgdGMF_reducedDims,
                                                        exprs_values=1,
                                                        assay.type=exprs_values)
{
  if(is.null(sgdGMF_reducedDims)){
    stop("first run 'runSGD' to obtain estimations for Imputation")
  }
  if(nrow(attr(sgdGMF_reducedDims,'rotation')) != nrow(x)){
      stop("all features should be used when performing 'runSGD' if imputation is wanted.")
  }
  imputedAssay <-.imputeMissingValues(assay(x, assay.type), sgdGMF_reducedDims, ...)
  assay(x, paste0(assay.type, '_imputed')) <- imputedAssay
  x
})

#' @export
#' @rdname SGDImpute
#' @importFrom SummarizedExperiment assay
setMethod("SGDImpute", "SingleCellExperiment", function(x,
                                                        reducedDimName = "SGD",
                                                        sgdGMF_reducedDims = reducedDim(x, reducedDimName),
                                                        ...,
                                                        exprs_values=1,
                                                        assay.type=exprs_values)
{
  if(is.null(sgdGMF_reducedDims)){
    stop("first run 'runSGD' to obtain estimations for Imputation")
  }
  if(nrow(attr(sgdGMF_reducedDims,'rotation')) != nrow(x)){
        stop("all features should be used when performing 'runSGD' if imputation is wanted.")
  }


  imputedAssay <- .imputeMissingValues(assay(x, assay.type), sgdGMF_reducedDims, ...)
  assay(x, paste0(assay.type, '_imputed')) <- imputedAssay
  x

})
