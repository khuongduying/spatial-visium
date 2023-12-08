library(SpatialExperiment)
library(spatialDE)
library(ggplot2)

# LOAD DATA
# Expression file used in python SpatialDE. 
data("Rep11_MOB_0")

# Sample Info file used in python SpatialDE
data("MOB_sample_info")

Rep11_MOB_0 <- Rep11_MOB_0[rowSums(Rep11_MOB_0) >= 3, ]

Rep11_MOB_0 <- Rep11_MOB_0[, row.names(MOB_sample_info)]
MOB_sample_info$total_counts <- colSums(Rep11_MOB_0)

head(MOB_sample_info)

X <- MOB_sample_info[, c("x", "y")]



# STABLILIZE
norm_expr <- stabilize(Rep11_MOB_0)

norm_expr[1:5, 1:5]

# regress_out
resid_expr <- regress_out(norm_expr, sample_info = MOB_sample_info)
resid_expr[1:5, 1:5]

# RUN
# For this example, run spatialDE on the first 1000 genes
sample_resid_expr <- head(resid_expr, 1000)

results <- spatialDE::run(sample_resid_expr, coordinates = X)
head(results[order(results$qval), ])


# MODEL SEARCH
de_results <- results[results$qval < 0.05, ]

ms_results <- model_search(
    sample_resid_expr,
    coordinates = X,
    de_results = de_results
)

# To show ms_results sorted on qvalue, uncomment the following line
head(ms_results[order(ms_results$qval), ])

head(ms_results)

# SPATIAL PATTERN
sp <- spatial_patterns(
    sample_resid_expr,
    coordinates = X,
    de_results = de_results,
    n_patterns = 4L, length = 1.5
)
sp$pattern_results

# VISUALIZE
gene <- "Rbfox3"

ggplot(data = MOB_sample_info, aes(x = x, y = y, color = norm_expr[gene, ])) +
    geom_point(size = 7) +
    ggtitle(gene) +
    scale_color_viridis_c() +
    labs(color = gene)


# Plot Spatial Patterns of Multiple Genes
ordered_de_results <- de_results[order(de_results$qval), ]

multiGenePlots(norm_expr,
    coordinates = X,
    ordered_de_results[1:6, ]$g,
    point_size = 4,
    viridis_option = "D",
    dark_theme = TRUE
)

FSV_sig(results, ms_results)
#> Warning: ggrepel: 8 unlabeled data points (too many overlaps). Consider
#> increasing max.overlaps


# SpatialExperiment integration

partial_counts <- head(Rep11_MOB_0, 1000)

spe <- SpatialExperiment(
    assays = list(counts = partial_counts),
    spatialData = DataFrame(MOB_sample_info[, c("x", "y")]),
    spatialCoordsNames = c("x", "y")
)

out <- spatialDE(spe, assay_type = "counts", verbose = FALSE)
head(out[order(out$qval), ])


# Plot Spatial Patterns of Multiple Genes (using SpatialExperiment object)

spe_results <- out[out$qval < 0.05, ]

ordered_spe_results <- spe_results[order(spe_results$qval), ]

multiGenePlots(spe,
    assay_type = "counts",
    ordered_spe_results[1:6, ]$g,
    point_size = 3,
    viridis_option = "D",
    dark_theme = TRUE
)

