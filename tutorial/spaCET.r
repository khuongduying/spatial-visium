library(SpaCET)

# PATH
fig_path <- "/mnt/rdisk/duydao/SPATIAL/GeneSmart/figures/tutorial/"
rdata <- "/mnt/rdisk/duydao/SPATIAL/GeneSmart/rdata/tutorial/"

# set the path to the in-house breast cancer ST data. User can set the paths to their own data.
visiumPath <- file.path(system.file(package = "SpaCET"), "extdata/Visium_BC")

# Create SpaCET object
## load ST data to create an SpaCET object.
SpaCET_obj <- create.SpaCET.object.10X(visiumPath = visiumPath)

## show this object.
str(SpaCET_obj)

# Show key quality control metrics
## calculate the QC metrics
SpaCET_obj <- SpaCET.quality.control(SpaCET_obj)

# plot the QC metrics
p1 <- SpaCET.visualize.spatialFeature(
    SpaCET_obj, 
    spatialType     = "QualityControl", 
    spatialFeatures = c("UMI","Gene"),
    imageBg         = TRUE
  )

ggsave(
    paste0(fig_path,"breast_QCmetrics.png"), 
    p1, 
    width   = 3000, 
    height  = 1500, 
    unit    = "px", 
    dpi     = 300
    )


# Deconvolve ST data
SpaCET_obj <- SpaCET.deconvolution(
                SpaCET_obj, 
                cancerType  = "BRCA", 
                coreNo      = 8
                )

saveRDS(
        SpaCET_obj,
        file          = paste0(rdata, "breast_spaCET_obj.rds"), 
        ascii         = FALSE, 
        version       = NULL,
        compress      = TRUE, 
        refhook       = NULL
        )

### import spaCET object
SpaCET_obj <- readRDS(paste0(rdata,"breast_spaCET_obj.rds"))

## show the ST deconvolution results
SpaCET_obj@results$deconvolution$propMat[1:13,1:6]


# Visualize the cell type proportion
## Show the spatial distribution of malignant cells and macrophages.
p2 <- SpaCET.visualize.spatialFeature(
      SpaCET_obj, 
      spatialType       = "CellFraction", 
      spatialFeatures   = c("Malignant","Macrophage")
)

ggsave(paste0(fig_path,"breast_cellType.png"), p2, width=3500, height=1500, unit="px", dpi=300)

## show the spatial distribution of all cell types.
p3 <- SpaCET.visualize.spatialFeature(
          SpaCET_obj, 
          spatialType           = "CellFraction", 
          spatialFeatures       = "All", 
          sameScaleForFraction  = TRUE,
          pointSize             = 0.1, 
          nrow                  = 5
          )

ggsave(paste0(fig_path,"breast_AllCellType.png"), p3, width=3500, height=3000, unit="px", dpi=300)


# Interactive visualization panel
SpaCET.visualize.spatialFeature(SpaCET_obj,interactive=TRUE)





####################################
# Estimate cell-cell interactions ##
####################################

# Find co-localized cell-type pairs

## calculate the cell-cell colocalization.
SpaCET_obj <- SpaCET.CCI.colocalization(SpaCET_obj)

## visualize the cell-cell colocalization.
p4 <- SpaCET.visualize.colocalization(SpaCET_obj)
ggsave(paste0(fig_path,"colocalization.png"), p4, width=4500, height=2000, unit="px", dpi=300)


# Analyze the L-R network enrichment within ST spots

## calculate the L-R network score across ST spots.
SpaCET_obj <- SpaCET.CCI.LRNetworkScore(SpaCET_obj,coreNo=32)

## visualize the L-R network score.
p5 <- SpaCET.visualize.spatialFeature(
                      SpaCET_obj, 
                      spatialType     = "LRNetworkScore", 
                      spatialFeatures = c("Network_Score","Network_Score_pv")
                    )
ggsave(paste0(fig_path,"LRNetwork.png"), p5, width=4500, height=2000, unit="px", dpi=300)


## Ligand-Receptor analysis for a co-localized cell-type pair
SpaCET_obj <- SpaCET.CCI.cellTypePair(
                          SpaCET_obj, 
                          cellTypePair = c("CAF","Macrophage M2")
                          )
### [1] "CAF and Macrophage M2 have potential intercellular interaction in the current tissue."

## Visualize the interaction analysis of a co-localized cell-type pair.
p6 <- SpaCET.visualize.cellTypePair(
                          SpaCET_obj, 
                          cellTypePair = c("CAF","Macrophage M2")
                          )

ggsave(paste0(fig_path,"interaction_colocalized_celltype.png"), p6, width=5500, height=2000, unit="px", dpi=300)
###


# Enrich cell-cell interactions at the tumor-immune interface
# source("/mnt/rdisk/duydao/SPATIAL/GeneSmart/scripts/bin/interaction_duy.r")
## Identify the Tumor-Stroma Interface
SpaCET_obj <- SpaCET.identify.interface(SpaCET_obj)

## Visualize the Interface
p7 <- SpaCET.visualize.spatialFeature(
                      SpaCET_obj, 
                      spatialType     = "Interface", 
                      spatialFeature  = "Interface"
                      )
ggsave(paste0(fig_path,"interface.png"), p7, width=3000, height=3000, unit="px", dpi=300)

## Combine the interface and interaction spots
## Import source code
# source("/mnt/rdisk/duydao/SPATIAL/GeneSmart/scripts/bin/interaction_duy.r")

SpaCET_obj <- SpaCET.combine.interface(SpaCET_obj, cellTypePair=c("CAF","Macrophage M2"))

## Visualize the Interface. The spatialFeature should be formated as "Interface&celltype1_celltype2". Celltype 1 and 2 is in the order of alphabet.
p8 <- SpaCET.visualize.spatialFeature(
                      SpaCET_obj, 
                      spatialType = "Interface", 
                      spatialFeature = "Interface&CAF_Macrophage M2"
                      )
ggsave(paste0(fig_path,"interface_celltype.png"), p8, width=3000, height=3000, unit="px", dpi=300)

## Compute the distance of CAF-M2 to tumor border
p9 <- SpaCET.distance.to.interface(
                      SpaCET_obj, 
                      cellTypePair = c("CAF", "Macrophage M2")
                      )
ggsave(paste0(fig_path,"distance_CAF_M2_tumorBorder.png"), p9, width=3000, height=3000, unit="px", dpi=300)



################################
## Explore cancer cell states ##
################################

## further deconvolve malignant cell states
SpaCET_obj <- SpaCET.deconvolution.malignant(SpaCET_obj, coreNo = 16)

## show cancer cell state fraction of the first five spots
SpaCET_obj@results$deconvolution$propMat[c("Malignant cell state A","Malignant cell state B"),1:6]

## Visualize 
p10 <- SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "CellFraction", 
  spatialFeatures=c("Malignant","Malignant cell state A","Malignant cell state B"), 
  nrow=1
)
ggsave(paste0(fig_path,"cancer_stages.png"), p10, width=6000, height=2000, unit="px", dpi=300)
