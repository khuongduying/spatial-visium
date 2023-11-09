library(BayesSpace)
library(ggplot2)

# PROCESSING THE DATA
melanoma <- getRDS("2018_thrane_melanoma", "ST_mel1_rep2")

set.seed(100)
dec <- scran::modelGeneVar(melanoma)
top <- scran::getTopHVGs(dec, n = 2000)

set.seed(101)
melanoma <- scater::runPCA(melanoma, subset_row = top)

## Add BayesSpace metadata
melanoma <- spatialPreprocess(melanoma, platform="Visium", skip.PCA=TRUE)

q <- 4  # Number of clusters
d <- 7  # Number of PCs

# CLUSTERING WITH BAYESSPACE
## Run BayesSpace clustering
set.seed(100)
melanoma <- spatialCluster(melanoma, q=q, d=d, platform='Visium',
                        nrep=50000, gamma=2)

## View results
palette <- c("purple", "red", "blue", "yellow", "darkblue")
clusterPlot(melanoma, palette=palette, color="black", size=0.1) +
    labs(title="BayesSpace")

# ENHANCING RESOLUTION WITH BAYESSPACE
set.seed(100)
melanoma.enhanced <- spatialEnhance(melanoma, q=q, d=d, platform="ST",
                                    nrep=50000, gamma=2, 
                                    verbose=TRUE, save.chain=TRUE,
                                    jitter_scale=3.5, jitter_prior=0.3)

clusterPlot(melanoma.enhanced, palette=palette, color="black", size=0.1) +
            labs(title="Enhanced clustering")


## Enhancement of marker gene expression
markers <- list()
markers[["Tumor"]] <- c("PMEL")
markers[["Fibroblast"]] <- c("COL1A1")
markers[["Macrophage"]] <- c("CD14", "FCGR1A", "FCGR1B")
markers[["B-cell"]] <- c("CD19", "MS4A1")
markers[["T-cell"]] <- c("CD2", "CD3D", "CD3E", "CD3G", "CD7")

markers

melanoma.enhanced <- enhanceFeatures(melanoma.enhanced, melanoma,
                                    model="xgboost",
                                    feature_names=purrr::reduce(markers, c),
                                    nrounds=0)


## Aggregated the expression of marker genes within each cell type
sum_counts <- function(sce, features) { ##features are markers
    if (length(features) > 1) {
    colSums(logcounts(sce)[features, ])
    } else {
    logcounts(sce)[features, ]
    }
}

spot_expr <- purrr::map(markers, function(xs) sum_counts(melanoma, xs))
enhanced_expr <- purrr::map(markers, function(xs) sum_counts(melanoma.enhanced, xs))


## Plot the spatial expression of each cell type’s markers 
## at spot-level and enhanced subspot resolution
library(patchwork)

plot_expression <- function(sce, expr, title) {
    featurePlot(sce, expr, color=NA) +
    viridis::scale_fill_viridis(option="C") +
    labs(title=title, fill="Log-normalized\nexpression")
}

plot_expression_comparison <- function(cell_type) {
    spot.plot <- plot_expression(melanoma, 
                                spot_expr[[cell_type]], 
                                "Spot")
    enhanced.plot <- plot_expression(melanoma.enhanced,
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


# DIFFERENTIAL EXPRESSION ANALYSIS OF SPATIAL CLUSTERS

## Using the same 2,000 HVGs previously computed for PCA, 
## excluding ribosomal
hvgs <- top[grep("^RP[LS]", top, invert=TRUE)]

melanoma.enhanced <- enhanceFeatures(melanoma.enhanced, melanoma, 
                                    model="xgboost",
                                    feature_names=hvgs,
                                    nrounds=0)

library(dplyr)
library(Seurat)
## Convert SCE to Seurat object and use BayesSpace cluster as identifier
sobj <- Seurat::CreateSeuratObject(counts=logcounts(melanoma.enhanced),
                                assay='Spatial',
                                meta.data=as.data.frame(colData(melanoma.enhanced)))

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