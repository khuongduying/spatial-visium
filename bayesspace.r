
# IMPORT
library(SingleCellExperiment)
library(SpatialExperiment)
library(ggplot2)
library(BayesSpace)

# LOADING DATA
melanoma <- getRDS(dataset="2018_thrane_melanoma", sample="ST_mel1_rep2")

# PRE-PROCESS DATA
set.seed(102)
melanoma <- spatialPreprocess(melanoma, platform="ST", 
                            n.PCs=7, n.HVGs=2000, log.normalize=FALSE)

# CLUSTERING

## SELECTING THE NUMBER OF CLUSTERS
melanoma <- qTune(melanoma, qs=seq(2, 10), platform="ST", d=7)
qPlot(melanoma)

## CLUSTERING WITH BAYESSPACE
set.seed(149)
melanoma <- spatialCluster(melanoma, q=4, platform="ST", d=7,
                        init.method="mclust", model="t", gamma=2,
                        nrep=1000, burn.in=100,
                        save.chain=TRUE)

## VISUALIZING SPATIAL CLUSTERS
clusterPlot(melanoma)

clusterPlot(melanoma, palette=c("purple", "red", "blue", "yellow"), color="black") +
            theme_bw() +
            xlab("Column") +
            ylab("Row") +
            labs(fill="BayesSpace\ncluster", title="Spatial clustering of ST_mel1_rep2")

# ENHANCED RESOLUTION
## Clustering at enhanced resolution
melanoma.enhanced <- spatialEnhance(melanoma, q=4, platform="ST", d=7,
                                    model="t", gamma=2,
                                    jitter_prior=0.3, jitter_scale=3.5,
                                    nrep=1000, burn.in=100,
                                    save.chain=TRUE)

## Plot enhaced resolution
clusterPlot(melanoma.enhanced)


## Enhancing the resolution of gene expression
markers <- c("PMEL", "CD2", "CD19", "COL1A1")
melanoma.enhanced <- enhanceFeatures(melanoma.enhanced, melanoma,
                                    feature_names=markers,
                                    nrounds=0)

logcounts(melanoma.enhanced)[markers, 1:5]

rowData(melanoma.enhanced)[markers, ]


## Visualizing enhanced gene expression
featurePlot(melanoma.enhanced, "PMEL")

### compare the spatial expression of the imputed marker genes
enhanced.plots <- purrr::map(markers, function(x) featurePlot(melanoma.enhanced, x))
patchwork::wrap_plots(enhanced.plots, ncol=2)

### compare to the spot-level expression
spot.plots <- purrr::map(markers, function(x) featurePlot(melanoma, x))
patchwork::wrap_plots(c(enhanced.plots, spot.plots), ncol=4)
