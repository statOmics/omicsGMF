
# This code is adapted from the scater package.

#' Wrapper functions to create plots for specific types of reduced dimension
#' results in a SingleCellExperiment object, similar as the \code{scater}
#' package.
#'
#' @param object A \link[SingleCellExperiment]{SingleCellExperiment-class}
#' object.
#' @param ... Additional arguments to pass to
#' \code{\link[scater]{plotReducedDim}}
#' from the \code{scater} package.
#' @param ncomponents Numeric scalar indicating the number of dimensions
#' components to (calculate and) plot. This can also be a numeric vector,
#' see \code{\link[scater]{plotReducedDim}} for details
#' @param dimred A string or integer scalar indicating the reduced dimension
#' result in \code{\link[SingleCellExperiment]{reducedDims}}(object) to plot.
#'
#' @details This is a wrapper around \code{\link[scater]{plotReducedDim}} that
#' uses the "GMF" slot from the \code{\link[SingleCellExperiment]{reducedDims}}
#' to obtain a dimensionality reduction plot.
#'
#' @return
#' A \code{\link[ggplot2]{ggplot}} object.
#' @name plotGMF
#' @seealso
#' \code{\link[scater]{plotReducedDim}}, for the underlying calculations.
#' \code{\link[scater]{plotPCA}}, for a similar wrapper.
#' @author Alexandre Segers
#'
#' @examples
#' example_sce <- mockSCE(ncells = 200, ngenes = 100)
#' example_sce <- runCVGMF(example_sce,
#'                          exprs_values="counts",
#'                          family = poisson(),
#'                          ncomponents = c(1:5))
#' example_sce <- runGMF(example_sce,
#'                       exprs_values="counts",
#'                       family = poisson(),
#'                       ncomponents = 3)
#' plotGMF(example_sce)
NULL

#' @export
#' @rdname plotGMF
#' @importFrom scater plotReducedDim
#' @import ggplot2
plotGMF <- function(object, ..., ncomponents=2, dimred = "GMF")
{
    scater::plotReducedDim(object, ncomponents = ncomponents, dimred = dimred, ...)
}
