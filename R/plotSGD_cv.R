#' Functions to create a scree plot for model selection.
#'
#' @param x Output of \link{sgdgmf.cv}, \link{calculateSGD_cv} or
#' \link{runSGD_cv}.
#' @param ... For the \code{plot_cv} generic, additional arguments to
#' pass to specific methods.
#' @param name String specifying the name to be used to obtain the
#'  cross-validation table object in the \code{\link{metadata}}. It is possible
#'  to specify multiple names if multiple cross-validations have been ran. 
#' @param method Function for summarization of the cross-validation over 
#' multiple folds. Default is mean. 
#' @param criteria The model selection criteria that is plotted. Default is
#' 'dev' (deviance residuals), but 'mae', 'mse', 'aic' and 'bic are possible.
#' @details This function plots a screeplot based on the output of
#' \link{runSGD_cv}, \link{calculateSGD_cv} or \link{sgdgmf.cv}.
#'
#' @return
#' A \link{plot} object.
#' @name plotSGD_cv
#' @seealso
#' \code{\link{runSGD_cv}}, to perform the cross-validation.
#' @author Alexandre Segers
#'
#' @examples
#' example_sce <- mockSCE(ncells = 200, ngenes = 100)
#' example_sce <- runSGD_cv(example_sce,
#'                          exprs_values="counts",
#'                          family = poisson(),
#'                          ncomponents = c(1:5))
#' plotSGD_cv(example_sce)
NULL


#' @importFrom stats aggregate
#' @import ggplot2
.plot_cv <- function (x, 
                      method = "mean",
                      criteria = "dev")
{
    if(is.null(x$method)){
        x$method <- "cv_SGD"
    }
    cv_summarized <- aggregate(x = x[,criteria], 
                               by = list("Method" = x$method, "ncomp" = x$ncomp), 
                               FUN = method) 
    
    ggplot(data = cv_summarized, aes(x = ncomp, y = x, col = Method, group = Method)) + 
        geom_point() + 
        geom_line() +
        theme_bw() +
        ylab(paste0("Out of sample ", criteria)) + 
        xlab("number of components")
    
}


#' @export
#' @rdname plotSGD_cv
setMethod("plot_cv", "ANY", .plot_cv)


#' @export
#' @rdname plotSGD_cv
#' @importFrom S4Vectors metadata
plotSGD_cv <-  function(x, ..., name = "cv_SGD")
{
    cv_list <- lapply(name, FUN = function(i){
        cbind(S4Vectors::metadata(x)[[i]], "method" = i)
    })
    
    cv_table <- do.call(rbind, cv_list)
    
    plot_cv(x = cv_table,
           ...)
}




