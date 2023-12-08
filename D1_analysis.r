

# LOAD SPE OBJECT
## IMPORT
library(SpatialExperiment)

## LOAD DATASET
sampleDir <- "/mnt/d4t/SPATIAL/GeneSmart/output/D1/test_output/HTI_spaceranger_count_D1/outs"
sampleNames <- "HTI_D1"
spe <- read10xVisium(samples = sampleDir, 
                        sample_id = sampleNames, 
                        type = "sparse", 
                        data = "filtered", 
                        images = "lowres", 
                        load = TRUE)

# QUALITY CONTROL

## OVERVIEW

### IMPORT
library(ggspavis)
fig_path <- "/mnt/d4t/SPATIAL/GeneSmart/figures/"

### PLOT SPATIAL COORDINATES (SPOTS)
p <- plotSpots(spe, size=0.1)
ggsave(paste(fig_path,"SpatialCoords.png"), p, width = 1000, height = 1000, units = "px", dpi=300)


## CALCULATE QC METRICS
### IMPORT
library("scater")

### SUBSET TO KEEP ONLY SPOTS OVER-TISSUE
spe <- spe[, colData(spe)$in_tissue == 1]
dim(spe)

### IDENTIFY MITOCHONDRIAL GENES
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$symbol)
table(is_mito)

rowData(spe)$symbol[is_mito] # view extracted mitochondrial gene names

### CALCULATE PER-SPOT QC METRICS AND STORE IN COLDATA
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
head(colData(spe))


## SELECTING THRESHOLDS

### LIBRARY SIZE

#### histogram of library sizes
png(paste(fig_path,"LibSizes.png"), width=800, height=800)
hist(colData(spe)$sum, breaks = 50)
dev.off()

#### select QC threshold for library size
#### check the number of spots below this threshold
qc_lib_size <- colData(spe)$sum < 600 ## Filter UMI<600 
table(qc_lib_size)
#### Add to colData
colData(spe)$qc_lib_size <- qc_lib_size

#### check spatial pattern of discarded spots
p <- plotQC(spe, type = "spots", 
            discard = "qc_lib_size")
ggsave(paste(fig_path,"qc_LibSize_spots.png"), p, width=1000, height=1000, unit="px", dpi=300)


### NUMBER OF EXPRESSED FEATURES

#### histogram of numbers of expressed genes
png(paste(fig_path, "n_expressed_genes.png"), height=1000, width=1000)
hist(colData(spe)$detected, breaks = 20)
dev.off()

#### select QC threshold for number of expressed genes
#### select a threshold of 400 expressed genes per spot
qc_detected <- colData(spe)$detected < 400
table(qc_detected)

#### check spatial pattern of discarded spots
p <- plotQC(spe, type = "spots", 
        discard = "qc_detected")
ggsave(paste(fig_path,"qc_GenesDetected_spots.png"), p, width=1000, height=1000, unit="px", dpi=300)



### PROPORTION OF MITOCHONDRIAL READS

#### histogram of mitochondrial read proportions
png(paste(fig_path, "hist_MTgenes.png"), height=1000, width=1000)
hist(colData(spe)$subsets_mito_percent, breaks = 50)
dev.off()

#### select QC threshold for mitochondrial read proportion
#### a threshold of 20% for the mitochondrial read proportion
qc_mito <- colData(spe)$subsets_mito_percent > 20
table(qc_mito)

#### Add qc_mito
colData(spe)$qc_mito <- qc_mito

#### check spatial pattern of discarded spots
p <- plotQC(spe, type = "spots", 
        discard = "qc_mito")
ggsave(paste(fig_path,"qc_mito.png"), p, width=1000, height=1000, unit="px", dpi=300)


### REMOVE LOW-QUALITY SPOTS

#### number of discarded spots for each metric
apply(cbind(qc_lib_size, qc_detected, qc_mito), 2, sum)

## qc_lib_size qc_detected     qc_mito 
##          2           1           0 

#### combined set of discarded spots
discard <- qc_lib_size | qc_detected | qc_mito
table(discard)

#### store in object
colData(spe)$discard <- discard

#### check spatial pattern of combined set of discarded spots
p <- plotQC(spe, type = "spots", 
        discard = "discard")
ggsave(paste(fig_path,"qc_discard.png"), p, width=1000, height=1000, unit="px", dpi=300)

#### remove combined set of low-quality spots
spe <- spe[, !colData(spe)$discard]
dim(spe)



# NORMALIZATION
## IMPORT
library(scran)

## CALCULATE LIBRARY SIZE FACTORS
spe <- computeLibraryFactors(spe)

summary(sizeFactors(spe))

hist(sizeFactors(spe), breaks=50)

## CALCULATE LOGCOUNTS AND STORE IN OBJECT
spe <- logNormCounts(spe)

## check
assayNames(spe)

dim(counts(spe))
dim(logcounts(spe))


# FEATURE SELECTION

## remove mitochondrial genes
spe <- spe[!is_mito, ]
dim(spe)

## fit mean-variance relationship
dec <- modelGeneVar(spe)

## visualize mean-variance relationship
fit <- metadata(dec)

png(paste(fig_path,"mean_variance.png"), width=1000, height=800)
plot(fit$mean, fit$var, 
        xlab = "mean of log-expression", 
        ylab = "variance of log-expression")
curve(fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)
dev.off()

## select top HVGs
top_hvgs <- getTopHVGs(dec, prop = 0.1)
length(top_hvgs)


##================================
## SPATIALLY VIRABLE GENES (SVGs)
##================================

# FIND SVGS WITH SPARK-X

## Identify top SVGs
## Plot top SVGs

###===========###
### SPATIALDE ###
###===========###
library("spatialDE")

## Transform counts
mx_counts <- as.matrix(counts(spe)) ## to matrix
df_counts <- data.frame(mx_counts) ## to dataFrame
df_counts <- tibble::rownames_to_column(df_counts, "geneID")

## Merge with gene_symbol
gene_symbol <- rowData(spe)
gene_symbol <- data.frame(gene_symbol)
gene_symbol <- tibble::rownames_to_column(gene_symbol,"geneID")

df_counts <- merge(gene_symbol,df_counts,
        #by.x="geneID",
        by.y="geneID",
        all.x=TRUE)

## Check duplicated
df_counts$symbol[which(duplicated(df_counts$symbol) == TRUE)]

## Some symbol are duplicated --> make them uniq by using "make.names()"
uniq_symbol <- make.names(df_counts$symbol, unique=TRUE)

## Add uniq symbol for df_counts rownames
rownames(df_counts) <- uniq_symbol

## Remove geneID and symbol column
df_counts <- df_counts[,-which(colnames(df_counts)==c("geneID","symbol"))]

## Return a matrix counts
g_counts <- as.matrix(df_counts)


## STABILIZE
norm_expr <- stabilize(g_counts)

## REGRESS_OUT

## Make library info
### lib_info - dataframe rownames=cell; x, y, total_counts

### Extract coord of cells
coords <- spatialCoords(spe)
colnames(coords) <- c("x","y") ## change column names
coords <- data.frame(coords) ## Turn to dataframe

### Extract count on cells
lib_size <- colData(spe)["total"]
colnames(lib_size) <- "total_counts"
lib_size <- data.frame(lib_size)

### Merge dataframe
sample_lib_info <- merge(coords,lib_size,
                by="row.names",
                all=TRUE)
## save file
fig_rdata <- "/mnt/d4t/SPATIAL/GeneSmart/rdata/"
write.csv(sample_lib_info, 
        file=paste(fig_rdata, "sample_lib_info.csv"),
        row.names=FALSE)


### regress_out calculation
resid_expr <- regress_out(norm_expr, sample_info = sample_lib_info)
resid_expr[1:5,1:5]

saveRDS()

# RUN SPATIALDE
de_results <- spatialDE::run(resid_expr, 
                coordinates=coords, 
                verbose = FALSE)

write.csv(de_results, 
        file=paste(fig_rdata,"de_results_raw.csv"),
        row.names=FALSE)

## Filter SVG that qval <0.05
de_results_filtered <- de_results[de_results$qval < 0.05, ]
write.csv(de_results_filtered, 
        file=paste(fig_rdata,"de_results_ft_qval005.csv"),
        row.names=FALSE)

# View the top 10 SVGs
View(head(de_results_filtered[order(de_results_filtered$qval), ],10))
top10_SVGs <- head(de_results_filtered[order(de_results_filtered$qval), ],10)
top50_SVGs <- head(de_results_filtered[order(de_results_filtered$qval), ],50)

write.csv(top50_SVGs, 
        file=paste(fig_rdata,"top50_SVGs.csv"),
        row.names=FALSE)
View(top50_SVGs)

## PLOT A GENE
gene <- "CDX2"

#p <- 
ggplot(data = sample_lib_info, aes(x = x, y = y, color = norm_expr[gene, ])) +
        geom_point(size = 5) +
        ggtitle(gene) +
        scale_color_viridis_c() +
        labs(color = gene)
ggsave(paste(fig_path,"PDL1.png"), p, width=2400, height=3000, unit="px", dpi=300)


## MODEL SEARCH
de_results_filtered <- de_results[de_results$qval < 0.05, ]

ms_results <- model_search(
        resid_expr,
        coordinates = coords,
        de_results = de_results_filtered
)

# To show ms_results sorted on qvalue, uncomment the following line
head(ms_results[order(ms_results$qval), ])

head(ms_results)

# SPATIAL PATTERN
sp <- spatial_patterns(
                resid_expr,
                coordinates = coords,
                de_results = de_results_filtered,
                n_patterns = 4L, length = 1.5
)
sp$pattern_results

## Plot Spatial Patterns of Multiple Genes
ordered_de_results <- de_results_filtered[order(de_results_filtered$qval),]

## TOP 20 SVGs
p <- multiGenePlots(norm_expr,
        coordinates = coords,
        ordered_de_results[17:20, ]$g,
        point_size = 1,
        viridis_option = "D",
        dark_theme = TRUE
)
ggsave(paste(fig_path,"top20_HVGs_5.png"), p, width=2400, height=3000, unit="px", dpi=300)

## T-CELL MARKERS
tcell_genes <- c("THEMIS","CD8B2","CD8A","CD3E","CD3D")
p <- multiGenePlots(norm_expr,
        coordinates = coords,
        tcell_genes,
        point_size = 1,
        viridis_option = "D",
        dark_theme = TRUE
)
ggsave(paste(fig_path,"tcell_markers.png"), p, width=2400, height=3000, unit="px", dpi=300)

## B-CELL MARKERS
bcell_genes <- c("MS4A1","CD19","CR2","CD79A","CD79B")
p <- multiGenePlots(norm_expr,
        coordinates = coords,
        bcell_genes,
        point_size = 1,
        viridis_option = "D",
        dark_theme = TRUE
)
ggsave(paste(fig_path,"bcell_markers.png"), p, width=2400, height=3000, unit="px", dpi=300)


## MACROPHAGES MARKERS
macrophage_genes <- c("MRC1","MSR1","CD163","FCN2","CD14")
p <- multiGenePlots(norm_expr,
        coordinates = coords,
        macrophage_genes,
        point_size = 1,
        viridis_option = "D",
        dark_theme = TRUE
)
ggsave(paste(fig_path,"macrophage_markers.png"), p, width=2400, height=3000, unit="px", dpi=300)


## NK CELL MARKERS
NKcell_genes <- c("KLRK1","KIR2DL4","XCL2","XCL1","SPON2")
p <- multiGenePlots(norm_expr,
        coordinates = coords,
        NKcell_genes,
        point_size = 1,
        viridis_option = "D",
        dark_theme = TRUE
)
ggsave(paste(fig_path,"NK_markers.png"), p, width=2400, height=3000, unit="px", dpi=300)

## Cancer subtype
cancer_type_genes <- c("KRT7","TTF1","IL12B","KRT5","KRT16","KRT17","TP63","NAPSA")
p <- multiGenePlots(norm_expr,
        coordinates = coords,
        cancer_type_genes,
        point_size = 1,
        viridis_option = "D",
        dark_theme = TRUE
)
ggsave(paste(fig_path,"cancer_type_markers.png"), p, width=2400, height=3000, unit="px", dpi=300)

## neuroendocrine
neuroen_genes <- c("NCAM1","SYP","CHGA")
p <- multiGenePlots(norm_expr,
        coordinates = coords,
        neuroen_genes,
        point_size = 1,
        viridis_option = "D",
        dark_theme = TRUE
)
ggsave(paste(fig_path,"markers_neuroendocrine.png"), p, width=2400, height=3000, unit="px", dpi=300)

## HER2 and MET
her2_met_genes <- c("ERBB2","MET")
p <- multiGenePlots(norm_expr,
        coordinates = coords,
        her2_met_genes,
        point_size = 2.5,
        viridis_option = "D",
        dark_theme = TRUE
)
ggsave(paste(fig_path,"markers_her2_met.png"), p, width=3000, height=3000, unit="px", dpi=300)

## Addition test
test_genes <- c("KRT20","CDX2","KLK3")
p <- multiGenePlots(norm_expr,
        coordinates = coords,
        test_genes,
        point_size = 2,
        viridis_option = "D",
        dark_theme = TRUE
)
ggsave(paste(fig_path,"markers_test.png"), p, width=2400, height=3000, unit="px", dpi=300)


## Plot Fraction Spatial Variance vs Q-value
p <- FSV_sig(de_results, ms_results)
ggsave(paste(fig_path,"FSV_plot.png"), p, width=2400, height=3000, unit="px", dpi=300)


#> Warning: ggrepel: 8 unlabeled data points (too many overlaps). Consider
#> increasing max.overlaps



##===============================
## DIMENSIONALITY REDUCTION #####
##===============================

## compute PCA
set.seed(123)
spe <- runPCA(spe, subset_row = top_hvgs)

reducedDimNames(spe)

## compute UMAP on top 50 PCs
set.seed(123)
spe <- runUMAP(spe, dimred = "PCA")

reducedDimNames(spe)

dim(reducedDim(spe, "UMAP"))

## update column names for easier plotting
colnames(reducedDim(spe, "UMAP"))
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)


## VISUALIZATIONS

library(ggspavis)

### plot top 2 PCA dimensions
plotDimRed(spe, type = "PCA")

### plot top 2 UMAP dimensions
plotDimRed(spe, type = "UMAP")


##============##
## CLUSTERING ##
##============##

#  Non-spatial clustering on HVGs

## graph-based clustering
set.seed(123)
k <- 10
g <- buildSNNGraph(spe, k = k, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
table(clus)



# plot clusters in spatial x-y coordinates
color <- c("#F0E442","#0072B2","#E69F00","#56B4E9","#D55E00","#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666")
plotSpots(spe, annotate="label", 
        palette=color,
        size=5)

# plot clusters in PCA reduced dimensions
plotDimRed(spe, type = "PCA", 
        annotate = "label", palette=color, size=3)

# plot clusters in UMAP reduced dimensions
plotDimRed(spe, type = "UMAP", 
        annotate = "label", palette = color, size=3)


## Spatially-aware clustering

##==============
## BAYESPACE ===
##==============

## IMPORT PACKAGE
library(BayesSpace)
library(SingleCellExperiment)
library(ggplot2)


## LOAD DATA

## PREPROCESSING DATA

### Update column name for rowData
colnames(rowData(spe)) <- "gene_name"

# Add gene_id to rowData
rowData(spe)$gene_id <- rownames(spe)

# check
rowData(spe)

# colData
colnames(colData(spe))[2] <- "row"
colnames(colData(spe))[3] <- "col"

# check
colData(spe)[1:5,]

## 
set.seed(100)
dec <- scran::modelGeneVar(spe)
top <- scran::getTopHVGs(dec, n = 2000)

top

set.seed(101)
spe <- scater::runPCA(spe, subset_row = top)

## Add BayesSpace metadata
spe <- spatialPreprocess(spe, platform="Visium", skip.PCA=TRUE)

## SELECTING THE NUMBER OF CLUSTERS
spe <- qTune(spe, qs=seq(2, 20), platform="Visium", d=15)
#p <- 

qPlot(spe)
ggsave(paste(fig_path,"qplot.png"), p, width = 1000, height = 1000, units = "px", dpi=300)

q <- 14  # Number of clusters
d <- 15  # Number of PCs
color <- c("#F0E442","#0072B2","#E69F00","#56B4E9","#C41DB0","#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#a9bc2d","#A6761D","#666666","#0000ff")


set.seed(100)
spe <- spatialCluster(spe, q=q, d=d, platform='Visium',
                        nrep=200000, gamma=3)

## View results
palette <- c("purple", "red", "blue", "yellow", "darkblue")
#p <- 

p <- clusterPlot(spe, palette=color, color="black", size=0.1) +
        labs(title="Cluster Plot of sample D1")
ggsave(paste(fig_path,"cluster_q14_i200K.png"), p, width = 2000, height = 2000, units = "px", dpi=300)

saveRDS(spe, file = paste(fig_rdata,"spatialD1.rds"), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

# ENHANCING RESOLUTION WITH BAYESSPACE
set.seed(100)
spe.enhanced <- spatialEnhance(spe, q=q, d=d, platform="Visium",
                                nrep=200000, gamma=3, 
                                verbose=TRUE, save.chain=TRUE,
                                jitter_scale=3.5, jitter_prior=0.3)

clusterPlot(spe.enhanced, palette=palette, color="black", size=0.1) +
        labs(title="Enhanced clustering")