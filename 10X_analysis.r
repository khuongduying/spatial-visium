library(SingleCellExperiment)
library(SpatialExperiment)
library(ggplot2)
library(BayesSpace)
library(scran)
library(scater)
library(ggspavis)
library(patchwork)

# PATH
fig_path <- "/mnt/d4t/SPATIAL/GeneSmart/figures/10X/"
rdata <- "/mnt/d4t/SPATIAL/GeneSmart/rdata/"


# LOAD DATA

## Output Spaceranger directory
sampleDir <- "/mnt/d4t/SPATIAL/GeneSmart/output/10X/Sample_spaceranger_count_10X/outs"

## CUSTOM READVISIUM
source("/mnt/d4t/SPATIAL/GeneSmart/scripts/readVisium_duy.r")

## Load to a Single-cell Experiment 
sce <- readVisium_duy(dirname = sampleDir)


# PRE-PROCESSING

## QUALITY CONTROL (QC)
library(scater)

### subset to keep only spots over tissue
sce <- sce[, colData(sce)$in_tissue == 1]

### identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(sce)$gene_name)

### calculate per-spot QC metrics
sce <- addPerCellQC(sce, subsets = list(mito = is_mito))

### select QC thresholds
qc_lib_size <- colData(sce)$sum < 600
qc_detected <- colData(sce)$detected < 400
qc_mito <- colData(sce)$subsets_mito_percent > 28
## qc_cell_count <- colData(sce)$cell_count > 10

### combined set of discarded spots
discard <- qc_lib_size | qc_detected | qc_mito
colData(sce)$discard <- discard

### filter low-quality spots
sce <- sce[, !colData(sce)$discard]



# NORMALIZATION

library(scran)
# calculate logcounts using library size factors
sce <- logNormCounts(sce)


# FEATURE SELECTION

## remove mitochondrial genes
sce <- sce[!is_mito, ]

## fit mean-variance relationship
dec <- modelGeneVar(sce)

## select top HVGs
top_hvgs <- getTopHVGs(dec, prop = 0.1)




############ TEST CLUSTERING SPOT LEVEL(SUCCESS)
## compute PCA for 2000 HVGs
set.seed(123)
sce.enhanced2 <- runPCA(sce.enhanced2, subset_row = top)

## compute UMAP on top 50 PCs
set.seed(123)
sce.enhanced2 <- runUMAP(sce.enhanced2, dimred = "PCA")

### plot top 2 PCA dimensions
plotDimRed(sce.enhanced2, type = "PCA")

# store cluster labels in column 'label' in colData
colLabels(sce.enhanced2) <- factor(colData(sce.enhanced2)$spatial.cluster)

### plot top 2 UMAP dimensions
p<-plotDimRed(sce.enhanced2, type = "UMAP",
          annotate = "label",
          palette = color, size=0.5) 
p
ggsave(file=paste(fig_path,"UMAP_q13_A1.png"),p, width=2000, height=3000, unit="px", dpi=300)
#################



# DIMENSIONALITY REDUCTION

## compute PCA
set.seed(123)
sce <- runPCA(sce, subset_row = top_hvgs)

## compute UMAP on top 50 PCs
set.seed(123)
sce <- runUMAP(sce, dimred = "PCA")

## update column names
colnames(reducedDim(sce, "UMAP")) <- paste0("UMAP", 1:2)

reducedDimNames(sce)


## VISUALIZATIONS
library(ggspavis)

### plot top 2 PCA dimensions
plotDimRed(sce, type = "PCA")

### plot top 2 UMAP dimensions
plotDimRed(sce, type = "UMAP") 


# BAYESSPACE

set.seed(100)
dec <- scran::modelGeneVar(sce)
top <- scran::getTopHVGs(dec, n = 2000)

set.seed(101)
sce <- scater::runPCA(sce, subset_row = top)

## Add BayesSpace metadata
sce <- spatialPreprocess(sce, platform="Visium", skip.PCA=TRUE)


## SELECTING THE NUMBER OF CLUSTERS
sce <- qTune(sce, qs=seq(2, 20), platform="Visium", d=15)

p <- qPlot(sce)
p
ggsave(file=paste(fig_path,"qPlot_10X.png"),p, width=2400, height=3000, unit="px", dpi=300)


q <- 15  # Number of clusters
d <- 15  # Number of PCs

# CLUSTERING WITH BAYESSPACE
## Run BayesSpace clustering
set.seed(100)
sce <- spatialCluster(sce, q=q, d=d, platform='Visium',
                        nrep=50000, gamma=3)

## View results
color <- c("#F0E442","#0072B2","#E69F00","#56B4E9","#D55E00","#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#e8affa","#ad0c22")
p <- clusterPlot(sce, palette=color, color="black", size=0.1) +
    labs(title="10X Clusters")
p
ggsave(paste(fig_path,"cluster_q15_d15_i50K.png"), p, width=3000, height=3000, unit="px", dpi=300)


# SAVE DATA
rdata <- "/mnt/d4t/SPATIAL/GeneSmart/rdata/"

saveRDS(sce, file = paste(rdata,"spatial10X_q15_ready.rds"), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)


# Import data
sce <- readRDS(paste(rdata,"spatialA1_q13_ready.rds"))

# ENHANCING RESOLUTION WITH BAYESSPACE
set.seed(100)
sce.enhanced <- spatialEnhance(sce, q=q, d=d, platform="Visium",
                                    nrep=50000, gamma=3, 
                                    verbose=TRUE, save.chain=TRUE,
                                    jitter_scale=3.5, jitter_prior=0.3)

## save sce.enhanced
saveRDS(sce.enhanced, file = paste(rdata,"spatial10X_q15_enhanced.rds"), ascii=FALSE, version=NULL,
        compress=TRUE, refhook=NULL)
## Import sce.enhanced
sce.enhanced <- readRDS("/mnt/d4t/SPATIAL/GeneSmart/rdata/ spatial10X_q15_enhanced.rds")


# Plot enhanced cluster

color <- c("#F0E442","#0072B2","#E69F00","#56B4E9","#D55E00","#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#e8affa","#ad0c22")


p <- clusterPlot(sce.enhanced, palette=color, color="lightgray", size=0.05) +
            labs(title="Enhanced clustering")
p
ggsave(paste(fig_path,"cluster_10X_enhanced_q15_i50K.png"), p, width=3000, height=3000, unit="px", dpi=300)


## Enhancement of marker gene expression
markers <- list()
markers[["Tumor"]] <- c("KRT7","TTF1","IL12B","KRT5","KRT16","KRT17","TP63","NAPSA")
markers[["Fibroblast"]] <- c("COL1A1")
markers[["Macrophage"]] <- c("MRC1","MSR1","CD163","FCN2","CD14")
markers[["B-cell"]] <- c("MS4A1","CD19","CR2","CD79A","CD79B")
markers[["T-cell"]] <- c("THEMIS","CD8B2","CD8A","CD3E","CD3D")
markers[["HER_MET"]] <- c("ERBB2","MET")
markers[["Neuroendocrine"]] <- c("NCAM1","SYP","CHGA")
markers[["NK-cell"]] <- c("KLRK1","KIR2DL4","XCL2","XCL1","SPON2")
markers

sce.enhanced <- enhanceFeatures(sce.enhanced, sce,
                                    model="xgboost",
                                    feature_names=purrr::reduce(markers, c),
                                    nrounds=0)


## Aggregated the expression of marker genes within each cell type
sum_counts <- function(SCE, features) { ##features are markers
    if (length(features) > 1) {
    colSums(logcounts(SCE)[features, ])
    } else {
    logcounts(SCE)[features, ]
    }
}

spot_expr <- purrr::map(markers, function(xs) sum_counts(sce, xs))
enhanced_expr <- purrr::map(markers, function(xs) sum_counts(sce.enhanced, xs))


## Plot the spatial expression of each cell type’s markers 
## at spot-level and enhanced subspot resolution
library(patchwork)

plot_expression <- function(SCE, expr, title) {
    featurePlot(SCE, expr, color=NA) +
    viridis::scale_fill_viridis(option="C") +
    labs(title=title, fill="Log-normalized\nexpression")
}

plot_expression_comparison <- function(cell_type) {
    spot.plot <- plot_expression(sce, 
                                spot_expr[[cell_type]], 
                                "Spot")
    enhanced.plot <- plot_expression(sce.enhanced,
                                    enhanced_expr[[cell_type]], 
                                    "Enhanced")
    (spot.plot + enhanced.plot) + 
    plot_annotation(title=cell_type,
                    theme=theme(plot.title=element_text(size=18)))
}
## Plot
plot_expression_comparison("Tumor")
plot_expression_comparison("Fibroblast")
plot_expression_comparison("Macrophage")
plot_expression_comparison("B-cell")
plot_expression_comparison("T-cell")

plot_expression_comparison("NK-cell")


p <- plot_expression_comparison("Tumor")
ggsave(paste(fig_path,"markers_tumor.png"), p, width=3000, height=1500, unit="px", dpi=300)

## Plot each marker gene


enhanced.plots <- purrr::map(
    markers[["NK-cell"]], function(x) featurePlot(sce.enhanced, x) + 
    viridis::scale_fill_viridis(option="C")
    )
p <- patchwork::wrap_plots(enhanced.plots, nrow=2, ncol=3)
p
ggsave(paste(fig_path,"markers_NKcell.png"), p, width=3000, height=1500, unit="px", dpi=300)


## Plot 1 gene
### Normal resolution
gene_name = "KRT17"
p<-featurePlot(sce, gene_name, color=NA) +
    viridis::scale_fill_viridis(option="C") +
    labs(title=gene_name, fill="Log-normalized\nexpression")
p
ggsave(paste(fig_path,"markers_KRT17_lowres.png"), p, width=3000, height=3000, unit="px", dpi=300)


### Enhanced resolution
gene_name = "KRT17"
p<-featurePlot(sce.enhanced, gene_name, color=NA) +
    viridis::scale_fill_viridis(option="C") +
    labs(title=gene_name, fill="Log-normalized\nexpression")
p
ggsave(paste(fig_path,"markers_KRT17_hires.png"), p, width=3000, height=3000, unit="px", dpi=300)


## SAVE DATA
### Path RDATA
rdata <- "/mnt/d4t/SPATIAL/GeneSmart/rdata/"

### spe
#saveRDS(spe, file = paste(rdata,"spatialD1_ready.rds"), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

### spe.enhanced
#saveRDS(spe.enhanced, file = paste(rdata,"spatialD1_enhanced_ready.rds"), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

### Check
hello <-readRDS(paste(rdata,"spatialD1_enhanced_ready.rds"))

# DIFFERENTIAL EXPRESSION ANALYSIS OF SPATIAL CLUSTERS

## Using the same 2,000 HVGs previously computed for PCA, 
## excluding ribosomal
hvgs <- top[grep("^RP[LS]", top, invert=TRUE)]

## enhaced top 2000 HVGs
sce.enhanced <- enhanceFeatures(sce.enhanced, sce, 
                                    model="xgboost",
                                    feature_names=hvgs,
                                    nrounds=0)

## save data
saveRDS(sce.enhanced, file = paste(rdata,"spatial10X_q15_enhanced_2000hvgs.rds"), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

library(dplyr)
library(Seurat)
## Convert SCE to Seurat object and use BayesSpace cluster as identifier
sobj <- Seurat::CreateSeuratObject(counts=logcounts(sce.enhanced),
                                assay='Spatial',
                                meta.data=as.data.frame(colData(sce.enhanced)))

sobj <- Seurat::SetIdent(sobj, value = "spatial.cluster")

## Scale data
sobj@assays$Spatial@scale.data <-
    sobj@assays$Spatial@data %>% as.matrix %>% t %>% scale %>% t

## Select top n markers from each cluster (by log fold change)
top_markers <- Seurat::FindAllMarkers(sobj, assay='Spatial', slot='data',
                                    group.by='spatial.cluster',
                                    only.pos=TRUE) %>% 
    group_by(cluster) %>% 
    top_n(5, avg_log2FC) ## problems

## Plot expression of markers
p<-Seurat::DoHeatmap(sobj, features = top_markers$gene, slot='scale.data',
                group.by = "spatial.cluster", group.colors=color, 
                angle=0, size=4, label = FALSE, raster=FALSE) + 
  guides(col = FALSE)
p
ggsave(paste(fig_path,"DEGs_q15.png"), p, width=2400, height=3000, unit="px", dpi=300)


# DIFFERENTIAL EXPRESSION OF TUMOR BORDER AND TUMOR-PROXIMAL LYMPHOID TISSUE
label_spot <- Vectorize(function(cluster, col) {
  if (cluster == 6) {
    "Cluster 6"
  } else if (cluster == 9) {
    "Cluster 9"
  } else {
    "Other"
  }
})

color <- c("#F0E442","#0072B2","#E69F00","#56B4E9","#D55E00","#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666")


DE.labels <- label_spot(sce.enhanced$spatial.cluster, sce.enhanced$col)
DE.labels <- factor(DE.labels, levels = c("Cluster 6", "Cluster 9", "Other"))
p<-clusterPlot(sce.enhanced, label=DE.labels, color="lightgray", size=0.1) + 
  scale_fill_manual(values=c('#1B9E77', '#E7298A', 'lightgray')) +
  labs(fill="Region")
p
ggsave(paste(fig_path,"cluster_6_vs_9.png"), p, width=3000, height=3000, unit="px", dpi=300)


## Perform DE analysis
## Add tumor proximal/border labels
sobj <- Seurat::AddMetaData(sobj, DE.labels, col.name="DE.label")
sobj <- Seurat::SetIdent(sobj, value = "DE.label")

## Subset to the two clusters of interest and re-scale
sobj <- subset(sobj, DE.label %in% c("Cluster 6", "Cluster 9"))
sobj@assays$Spatial@scale.data <-
  sobj@assays$Spatial@data %>% as.matrix %>% t %>% scale %>% t

## Select top n markers from each cluster (by log fold change)
top_markers <- Seurat::FindAllMarkers(sobj, assay='Spatial', slot='data',
                                      group.by='DE.label',
                                      only.pos=TRUE) %>% 
  group_by(cluster) %>% 
  top_n(5, avg_log2FC)

## Plot expression of markers
p<-Seurat::DoHeatmap(sobj, features = top_markers$gene, slot='scale.data',
                  group.by = "DE.label", group.colors=c('#1B9E77', '#E7298A'), 
                  angle=0, size=4, label = FALSE, raster=FALSE) + 
  guides(col = FALSE)
p
ggsave(paste(fig_path,"DEGs_cluster_6_vs_9.png"), p, width=3000, height=3000, unit="px", dpi=300)
