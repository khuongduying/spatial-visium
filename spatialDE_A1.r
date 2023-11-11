###===========###
### SPATIALDE ###
###===========###
library("spatialDE")

## LOAD DATASET
sampleDir <- "/mnt/d4t/SPATIAL/GeneSmart/output/A1/HTI_spaceranger_count_A1/outs"
sampleNames <- "HTI_A1"
spe <- read10xVisium(samples = sampleDir, 
                        sample_id = sampleNames, 
                        type = "sparse", 
                        data = "filtered", 
                        images = "lowres", 
                        load = TRUE)

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
        file=paste(fig_rdata, "A1_sample_lib_info.csv"),
        row.names=FALSE)


### regress_out calculation
resid_expr <- regress_out(norm_expr, sample_info = sample_lib_info)
resid_expr[1:5,1:5]

saveRDS(resid_expr, file = paste0(fig_rdata,"resid_expr_A1.rds"), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

saveRDS(coords, file = paste0(fig_rdata,"coords_A1.rds"), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

resid_expr <- readRDS(paste0(fig_rdata,"resid_expr_A1.rds"))
coords <- readRDS(paste0(fig_rdata,"coords_A1.rds"))
# RUN SPATIALDE
de_results <- spatialDE::run(resid_expr, 
                coordinates=coords, 
                verbose = FALSE)

write.csv(de_results, 
        file=paste(fig_rdata,"de_results_raw_A1.csv"),
        row.names=FALSE)

## Filter SVG that qval <0.05
de_results_filtered <- de_results[de_results$qval < 0.05, ]
write.csv(de_results_filtered, 
        file=paste(fig_rdata,"de_results_ft_qval005_A1.csv"),
        row.names=FALSE)

# View the top 10 SVGs
View(head(de_results_filtered[order(de_results_filtered$qval), ],10))
top10_SVGs <- head(de_results_filtered[order(de_results_filtered$qval), ],10)
top50_SVGs <- head(de_results_filtered[order(de_results_filtered$qval), ],50)

write.csv(top50_SVGs, 
        file=paste(fig_rdata,"top50_SVGs_A1.csv"),
        row.names=FALSE)
View(top50_SVGs)