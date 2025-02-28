#' Impute missing values based on the results of runSGD.
#'
#' @param x a numeric matrix of expression counts or mass spectrometry intensities containing missing values and with
#' features in the rows and samples in columns.
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such a matrix.
#' @param reducedDimName the name of the \code{\link{reducedDim}} slot corresponding to the dimensionality reduction
#' obtained with runSGD when \code{x} is a \linkS4class{SingleCellExperiment}.
#' @param sgdGMF_reducedDims the output obtained by \code{\link{runSGD}} or \code{\link{calculateSGD}}. If \code{x} is a
#' \linkS4class{SingleCellExperiment}, \code{sgdGMF_reducedDims} is taken from \code{\link{reducedDim}(x, 'SGD')}.
#' @param assay.type Integer scalar or string indicating which assay of \code{x} contains the values of interest.
#' @param exprs_values Alias to \code{assay.type}.
#' @param ... For the \code{SGDImpute} generic, additional arguments to pass to specific methods.
#'
#' @details
#' Imputation is only possible after running \code{runSGD} using all features.
#' Therefore, \code{subset_row} or \code{ntop} should be set to NULL.
#'
#' @return
#' A SingleCellExperiment object is returned with an extra assay containing the imputed values.
#'
#' @name SGDImpute
#' @seealso
#' \code{\link{runSGD}}, to conveniently obtain the matrix factorization.
#'
#' @author Alexandre Segers
#'
#' @examples
#' example_sce <- mockSCE(ncells = 200, ngenes = 100)
#' example_sce <- logNormCounts(example_sce)
#' assay(example_sce, 'logcounts')[assay(example_sce, 'logcounts') == 0] <- NA
#' example_sce <- runSGD(example_sce,
#'                       exprs_values="logcounts",
#'                       family = gaussian(),
#'                       ncomponents = 3)
#' example_sce <- SGDImpute(example_sce)
NULL

.imputeMissingValues <- function(x, sgdGMF_reducedDims)
{
    if(attr(sgdGMF_reducedDims,"family")$family != "gaussian"){
        stop("Imputation currently only possible for Gaussian family.")
    }

    if(is.null(attr(sgdGMF_reducedDims, "offset"))){
        means <- attr(sgdGMF_reducedDims,'B') %*% t(attr(sgdGMF_reducedDims,'X')) +
            attr(sgdGMF_reducedDims,'Z') %*% t(attr(sgdGMF_reducedDims,'Gamma')) +
            attr(sgdGMF_reducedDims, 'rotation') %*% t(sgdGMF_reducedDims)
    } else{
        means <- attr(sgdGMF_reducedDims,'B') %*% t(attr(sgdGMF_reducedDims,'X')) +
            attr(sgdGMF_reducedDims,'Z') %*% t(attr(sgdGMF_reducedDims,'Gamma')) +
            attr(sgdGMF_reducedDims, 'rotation') %*% t(sgdGMF_reducedDims) +
            attr(sgdGMF_reducedDims,'offset')
    }

    x[is.na(x)] <- means[is.na(x)]
    x
}


#' @export
#' @rdname SGDImpute
setMethod("SGDImpute", "ANY", .imputeMissingValues)




#' @export
#' @rdname SGDImpute
#' @importFrom SummarizedExperiment assay assay<-
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

  assay(x, paste0(assay.type, '_imputed')) <-  .imputeMissingValues(assay(x, assay.type), sgdGMF_reducedDims)
  x
})

#' @export
#' @rdname SGDImpute
#' @importFrom SummarizedExperiment assay assay<-
setMethod("SGDImpute", "SingleCellExperiment", function(x,
                                                        reducedDimName = "SGD",
                                                        sgdGMF_reducedDims = reducedDim(x, reducedDimName),
                                                        exprs_values=1,
                                                        assay.type=exprs_values)
{
  if(is.null(sgdGMF_reducedDims)){
    stop("first run 'runSGD' to obtain estimations for Imputation")
  }
  if(nrow(attr(sgdGMF_reducedDims, 'rotation')) != nrow(x)){
        stop("all features should be used when performing 'runSGD' if imputation is wanted.")
  }

  assay(x, paste0(assay.type, '_imputed')) <- .imputeMissingValues(assay(x, assay.type), sgdGMF_reducedDims)
  x

})
