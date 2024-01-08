library(SingleCellExperiment)
library(SpatialExperiment)
library(ggplot2)
library(BayesSpace)
library(scran)
library(scater)

## There are errors in the source of readVisium() function
sampleDir <- "/mnt/d4t/SPATIAL/GeneSmart/output/D1/test_output/HTI_spaceranger_count_D1/outs"
spe <- readVisium(dirname = sampleDir)

## Use this modify function to READ VISIUM DATA
source('/mnt/d4t/SPATIAL/GeneSmart/scripts/readVisium_duy.r')

# IMPORT DATA
spe <- readVisium_duy(sampleDir)

## Or read from the saved data
spe <- readRDS("/mnt/d4t/SPATIAL/GeneSmart/rdata/ spatialD1_ready.rds")
spe.enhanced <- readRDS("/mnt/d4t/SPATIAL/GeneSmart/rdata/ spatialD1_enhanced_ready.rds")

# BAYESSPACE
set.seed(100)
dec <- scran::modelGeneVar(spe)
top <- scran::getTopHVGs(dec, n = 2000)

set.seed(101)
spe <- scater::runPCA(spe, subset_row = top)

## Add BayesSpace metadata
spe <- spatialPreprocess(spe, platform="Visium", skip.PCA=TRUE)

q <- 4  # Number of clusters
d <- 15  # Number of PCs

# CLUSTERING WITH BAYESSPACE
## Run BayesSpace clustering
set.seed(100)
spe <- spatialCluster(spe, q=q, d=d, platform='Visium',
                        nrep=50000, gamma=3)

## View results
palette <- c("purple", "red", "blue", "yellow", "darkblue")
p <- clusterPlot(spe, palette=color, color="black", size=0.1) +
    labs(title="BayesSpace")
ggsave(paste(fig_path,"cluster_q4_i50K.png"), p, width=2400, height=3000, unit="px", dpi=300)

# ENHANCING RESOLUTION WITH BAYESSPACE
set.seed(100)
spe.enhanced2 <- spatialEnhance(spe, q=q, d=d, platform="Visium",
                                    nrep=50000, gamma=3, 
                                    verbose=TRUE, save.chain=TRUE,
                                    jitter_scale=3.5, jitter_prior=0.3)

color <- c("#F0E442","#0072B2","#E69F00","#56B4E9","#D55E00","#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666")

p <- clusterPlot(spe.enhanced, palette=color, color="lightgray", size=0.05) +
            labs(title="Enhanced clustering")
p
ggsave(paste(fig_path,"cluster_enhanced_q10_i50K.png"), p, width=2400, height=3000, unit="px", dpi=300)


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

spe.enhanced <- enhanceFeatures(spe.enhanced, spe,
                                    model="xgboost",
                                    feature_names=purrr::reduce(markers, c),
                                    nrounds=0)


## Aggregated the expression of marker genes within each cell type
sum_counts <- function(SPE, features) { ##features are markers
    if (length(features) > 1) {
    colSums(logcounts(SPE)[features, ])
    } else {
    logcounts(SPE)[features, ]
    }
}

spot_expr <- purrr::map(markers, function(xs) sum_counts(spe, xs))
enhanced_expr <- purrr::map(markers, function(xs) sum_counts(spe.enhanced, xs))


## Plot the spatial expression of each cell type’s markers 
## at spot-level and enhanced subspot resolution
library(patchwork)

plot_expression <- function(SPE, expr, title) {
    featurePlot(SPE, expr, color=NA) +
    viridis::scale_fill_viridis(option="C") +
    labs(title=title, fill="Log-normalized\nexpression")
}

plot_expression_comparison <- function(cell_type) {
    spot.plot <- plot_expression(spe, 
                                spot_expr[[cell_type]], 
                                "Spot")
    enhanced.plot <- plot_expression(spe.enhanced,
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
plot_expression_comparison("NK-cell")


p <- plot_expression_comparison("Tumor")
ggsave(paste(fig_path,"markers_before_after_enhanced_2.png"), p, width=3000, height=1500, unit="px", dpi=300)

## Plot each marker gene

enhanced.plots <- purrr::map(
    markers[["Neuroendocrine"]], function(x) featurePlot(spe.enhanced, x) + 
    viridis::scale_fill_viridis(option="C") # or K
    )
p <- patchwork::wrap_plots(enhanced.plots, nrow=1, ncol=3)
p
ggsave(paste(fig_path,"markers_neuroendo_enhanced.png"), p, width=3000, height=1500, unit="px", dpi=300)





## SAVE DATA
### Path RDATA
rdata <- "/mnt/d4t/SPATIAL/GeneSmart/rdata/"

### spe
saveRDS(spe, file = paste(rdata,"spatialD1_ready.rds"), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

### spe.enhanced
saveRDS(spe.enhanced, file = paste(rdata,"spatialD1_enhanced_ready.rds"), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

### Check
hello <-readRDS(paste(rdata,"spatialD1_enhanced_ready.rds"))

# DIFFERENTIAL EXPRESSION ANALYSIS OF SPATIAL CLUSTERS

## Using the same 2,000 HVGs previously computed for PCA, 
## excluding ribosomal
hvgs <- top[grep("^RP[LS]", top, invert=TRUE)]

spe.enhanced <- enhanceFeatures(spe.enhanced, spe, 
                                    model="xgboost",
                                    feature_names=hvgs,
                                    nrounds=0)

library(dplyr)
library(Seurat)
## Convert SCE to Seurat object and use BayesSpace cluster as identifier
sobj <- Seurat::CreateSeuratObject(counts=logcounts(spe.enhanced),
                                assay='Spatial',
                                meta.data=as.data.frame(colData(spe.enhanced)))

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
Seurat::DoHeatmap(sobj, features = top_markers$gene, slot='scale.data',
                group.by = "spatial.cluster", group.colors=palette, 
                angle=0, size=4, label = FALSE, raster=FALSE) + 
  guides(col = FALSE)


# DIFFERENTIAL EXPRESSION OF TUMOR BORDER AND TUMOR-PROXIMAL LYMPHOID TISSUE
label_spot <- Vectorize(function(cluster, col) {
  if (cluster == 1 && col < 19) {
    "Tumor border"
  } else if (cluster == 4 && col < 19) {
    "Lymphoid"
  } else {
    "Other"
  }
})

DE.labels <- label_spot(melanoma.enhanced$spatial.cluster, melanoma.enhanced$col)
DE.labels <- factor(DE.labels, levels = c("Lymphoid", "Tumor border", "Other"))
clusterPlot(melanoma.enhanced, label=DE.labels, color="black", size=0.1) + 
  scale_fill_manual(values=c('#0173b2', '#de8f05', '#949494')) +
  labs(fill="Region")

## Perform DE analysis
## Add tumor proximal/border labels
sobj <- Seurat::AddMetaData(sobj, DE.labels, col.name="DE.label")
sobj <- Seurat::SetIdent(sobj, value = "DE.label")

## Subset to the two clusters of interest and re-scale
sobj <- subset(sobj, DE.label %in% c("Lymphoid", "Tumor border"))
sobj@assays$Spatial@scale.data <-
  sobj@assays$Spatial@data %>% as.matrix %>% t %>% scale %>% t

## Select top n markers from each cluster (by log fold change)
top_markers <- Seurat::FindAllMarkers(sobj, assay='Spatial', slot='data',
                                      group.by='DE.label',
                                      only.pos=TRUE) %>% 
  group_by(cluster) %>% 
  top_n(5, avg_log2FC)

## Plot expression of markers
Seurat::DoHeatmap(sobj, features = top_markers$gene, slot='scale.data',
                  group.by = "DE.label", group.colors=c('#0173b2', '#de8f05'), 
                  angle=0, size=4, label = FALSE, raster=FALSE) + 
  guides(col = FALSE)

############

library(SpatialExperiment)
library(scater)


example(read10xVisium, echo = FALSE)
spe

# PROCESSING THE DATA
#melanoma <- getRDS("2018_thrane_melanoma", "ST_mel1_rep2")

rownames(spe) <- rowData(spe)$symbol
spe


## Log counts
spe <- logNormCounts(spe)

set.seed(100)
dec <- scran::modelGeneVar(spe)
top <- scran::getTopHVGs(dec, n = 2000)

set.seed(101)
melanoma <- scater::runPCA(melanoma, subset_row = top)

## Add BayesSpace metadata
melanoma <- spatialPreprocess(melanoma, platform="ST", skip.PCA=TRUE)

q <- 4  # Number of clusters
d <- 7  # Number of PCs