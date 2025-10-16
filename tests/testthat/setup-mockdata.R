# Creating a mock SingleCellExperiment
set.seed(100)

sce <- scuttle::mockSCE(ncells = 20, ngenes = 10)

# Removing alternative experiments
SingleCellExperiment::altExps(sce) <- NULL

# Create Gaussian experiment
sce <- scater::logNormCounts(sce)

# Create proteomics data with missing values
assay(sce, "logintensities") <- assay(sce, "logcounts")
assay(sce, "logintensities")[assay(sce, "logintensities") == 0] <- NA