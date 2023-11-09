dirname <- "/mnt/d4t/SPATIAL/GeneSmart/output/D1/test_output/HTI_spaceranger_count_D1/outs"

readVisium_duy <- function(dirname) {
    spatial_dir <- file.path(dirname, "spatial")
        matrix_dir <- file.path(dirname, "filtered_feature_bc_matrix")
        
        if (!dir.exists(matrix_dir))
            stop("Matrix directory does not exist:\n  ", matrix_dir)
        if (!dir.exists(spatial_dir))
            stop("Spatial directory does not exist:\n  ", spatial_dir)
        
        colData <- read.csv(file.path(spatial_dir, "tissue_positions_list.csv"), header=TRUE)
        
        ## We're using spatialLIBD conventions here (legacy), but should eventually
        ##   update to canonical spaceRanger names:
        ##   c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
        ##   https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/images
        colnames(colData) <- c("spot", "in_tissue", "row", "col", "imagerow", "imagecol")
        rownames(colData) <- colData$spot
        colData <- colData[colData$in_tissue > 0, ]
        
        rowData <- read.table(file.path(matrix_dir, "features.tsv.gz"), header=FALSE)
        colnames(rowData) <- c("gene_id", "gene_name", "feature_type")
        rowData <- rowData[, c("gene_id", "gene_name")]
        rownames(rowData) <- scater::uniquifyFeatureNames(rowData$gene_id, rowData$gene_name)
        
        counts <- Matrix::readMM(file.path(matrix_dir, "matrix.mtx.gz"))
        barcodes <- read.table(file.path(matrix_dir, "barcodes.tsv.gz"), header=FALSE)
        colnames(counts) <- barcodes$V1
        rownames(counts) <- rownames(rowData)
        counts <- counts[, rownames(colData)]
        
        sce <- SingleCellExperiment(assays=list(counts=counts),
                                    rowData=rowData,
                                    colData=colData)
        
        metadata(sce)$BayesSpace.data <- list()
        metadata(sce)$BayesSpace.data$platform <- "Visium"
        metadata(sce)$BayesSpace.data$is.enhanced <- FALSE

        sce
}