
test_that("calculateGMF() works on poisson data", {
    set.seed(100)
    res1 <- reducedDim(runGMF(sce, ncomponents = 5), "GMF")
    set.seed(100)
    res2 <- calculateGMF(sce, ncomponents = 5)
    set.seed(100)
    res3 <- calculateGMF(assay(sce, 'counts'), ncomponents = 5)


    expect_identical(c(res1), c(res2), c(res3))
    expect_identical(attr(res1, "X"), attr(res2, "X"), attr(res3, "X"))
    expect_identical(attr(res1, "Beta"), attr(res2, "Beta"), attr(res3, "Beta"))
    expect_identical(attr(res1, "Beta"), attr(res2, "Beta"), attr(res3, "Beta"))
    expect_identical(attr(res1, "Z"), attr(res2, "Z"), attr(res3, "Z"))
    expect_identical(attr(res1, "Gamma"), attr(res2, "Gamma"), attr(res3, "Gamma"))
    expect_identical(attr(res1, "offset"), attr(res2, "offset"), attr(res3, "offset"))

})


test_that("calculateGMF() works when introducing X and Z matrices", {
    set.seed(100)
    X = matrix(rnorm(ncol(sce)*3), ncol = 3)
    Z = matrix(rnorm(nrow(sce)*3), ncol = 3)
    set.seed(100)
    res1 <- reducedDim(runGMF(sce, X = X, Z = Z, ncomponents = 5), "GMF")
    set.seed(100)
    res2 <- calculateGMF(sce, X = X, Z = Z, ncomponents = 5)
    set.seed(100)
    res3 <- calculateGMF(assay(sce, 'counts'), X = X, Z = Z, ncomponents = 5)


    expect_identical(c(res1), c(res2), c(res3))
    expect_identical(attr(res1, "X"), attr(res2, "X"), attr(res3, "X"))
    expect_identical(attr(res1, "Beta"), attr(res2, "Beta"), attr(res3, "Beta"))
    expect_identical(attr(res1, "Beta"), attr(res2, "Beta"), attr(res3, "Beta"))
    expect_identical(attr(res1, "Z"), attr(res2, "Z"), attr(res3, "Z"))
    expect_identical(attr(res1, "Gamma"), attr(res2, "Gamma"), attr(res3, "Gamma"))
    expect_identical(attr(res1, "offset"), attr(res2, "offset"), attr(res3, "offset"))


    expect_equal(dim(attr(res1, "X"))[2], 3)
    expect_equal(dim(attr(res1, "Z"))[2], 3)

})

test_that("calculateGMF() works on gaussian data with missing values", {
    set.seed(100)
    res1 <- reducedDim(runGMF(sce,
                              exprs_values = "logintensities",
                              family = gaussian(),
                              ncomponents = 5), "GMF")
    set.seed(100)
    res2 <- calculateGMF(sce,
                         exprs_values = "logintensities",
                         family = gaussian(),
                         ncomponents = 5)
    set.seed(100)
    res3 <- calculateGMF(assay(sce, 'logintensities'),
                         family = gaussian(),
                         ncomponents = 5)


    expect_identical(c(res1), c(res2), c(res3))
    expect_identical(attr(res1, "X"), attr(res2, "X"), attr(res3, "X"))
    expect_identical(attr(res1, "Beta"), attr(res2, "Beta"), attr(res3, "Beta"))
    expect_identical(attr(res1, "Beta"), attr(res2, "Beta"), attr(res3, "Beta"))
    expect_identical(attr(res1, "Z"), attr(res2, "Z"), attr(res3, "Z"))
    expect_identical(attr(res1, "Gamma"), attr(res2, "Gamma"), attr(res3, "Gamma"))
    expect_identical(attr(res1, "offset"), attr(res2, "offset"), attr(res3, "offset"))

})

test_that("calculateCVGMF() works on poisson data", {
    set.seed(100)
    res1 <- metadata(runCVGMF(sce, ncomponents = c(1:3)))$cv_GMF
    set.seed(100)
    res2 <- calculateCVGMF(sce, ncomponents = c(1:3))
    expect_identical(res1, res2)
    set.seed(100)
    res3 <- calculateCVGMF(assay(sce, 'counts'), ncomponents = c(1:3))
    expect_identical(res2, res3)


})

test_that("calculateCVGMF() works on gaussian data with missing values", {
    set.seed(100)
    res1 <- metadata(runCVGMF(sce,
                                    exprs_values = "logintensities",
                                    family = gaussian(),
                                    ncomponents = c(1:3)))$cv_GMF
    set.seed(100)
    res2 <- calculateCVGMF(sce,
                      exprs_values = "logintensities",
                      family = gaussian(),
                      ncomponents = c(1:3))
    expect_identical(res1, res2)
    set.seed(100)
    res3 <- calculateCVGMF(assay(sce, 'logintensities'),
                                           family = gaussian(),
                                           ncomponents = c(1:3))
    expect_identical(res2, res3)
})




test_that("calculateRankGMF() works on poisson data", {
    set.seed(100)
    res1 <- metadata(runRankGMF(sce, maxcomp = 5))$rank_GMF
    set.seed(100)
    res2 <- calculateRankGMF(sce, maxcomp = 5)
    expect_identical(res1, res2)
    set.seed(100)
    res3 <- calculateRankGMF(assay(sce, 'counts'), maxcomp = 5)
    expect_identical(res2, res3)

})



