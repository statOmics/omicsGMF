
test_that("imputeGMF() works on gaussian data with missing values", {
    set.seed(100)
    dimred <- runGMF(sce, exprs_values = "logintensities",
                     family = gaussian(),
                     ncomponents = 5)
    res1 <- assay(imputeGMF(dimred, exprs_values = "logintensities"),
                  "imputedAssay")
    res2 <- assay(imputeGMF(sce,
                      sgdGMF_reducedDims = reducedDim(dimred, "GMF"),
                      exprs_values = "logintensities"),
                  "imputedAssay")
    res3 <- imputeGMF(assay(sce, 'logintensities'),
                      sgdGMF_reducedDims = reducedDim(dimred, "GMF"))

    expect_identical(res1, res2, res3)
    expect_equal(sum(is.na(res1)), 0)
})
