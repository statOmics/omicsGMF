test_that("plotGMF() works", {
    set.seed(100)
    res1 <- runGMF(sce, ncomponents = 5)
    expect_s3_class(plotGMF(res1), "ggplot")
})


