
library(CARD)
library(MuSiC)

## PATH
test_CARD <- "/mnt/d4t/SPATIAL/GeneSmart/datasets/test_CARD/"

#### load the example spatial transcriptomics count data, 
load(paste0(test_CARD,"spatial_count.RData"))

#### load the example spatial location data, 
load(paste0(test_CARD, "spatial_location.RData"))

#### load the example scRNA-seq count data, 
load(paste0(test_CARD, "sc_count.RData"))

#### load the example scRNA-seq meta data, 
load(paste0(test_CARD, "sc_meta.RData"))

#### load the Marker list
load(paste0(test_CARD, "markerList.RData"))


# CREATE AN CARD OBJECT
CARD_obj = createCARDObject(
	sc_count = sc_count,
	sc_meta = sc_meta,
	spatial_count = spatial_count,
	spatial_location = spatial_location,
	ct.varname = "cellType",
	ct.select = unique(sc_meta$cellType),
	sample.varname = "sampleInfo",
	minCountGene = 100,
	minCountSpot = 5) 

# DECONVOLUTION USING CARD
CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)

print(CARD_obj@Proportion_CARD[1:2,])

## VISUALIZE THE PROPORTION OF EACH CELL TYPE

### set the colors. Here, I just use the colors in the manuscript, if the color is not provided, the function will use default color in the package. 
colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
    "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
    "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")
p1 <- CARD.visualize.pie(
	proportion = CARD_obj@Proportion_CARD,
	spatial_location = CARD_obj@spatial_location, 
 	colors = colors, 
  	radius = 0.4) ### You can choose radius = NULL or your own radius number
print(p1)


### select the cell type that we are interested
ct.visualize = c("Acinar_cells","Cancer_clone_A","Cancer_clone_B","Ductal_terminal_ductal_like","Ductal_CRISP3_high-centroacinar_like","Ductal_MHC_Class_II","Ductal_APOL1_high-hypoxic","Fibroblasts")
### visualize the spatial distribution of the cell type proportion
p2 <- CARD.visualize.prop(
	proportion = CARD_obj@Proportion_CARD,        
	spatial_location = CARD_obj@spatial_location, 
	ct.visualize = ct.visualize,                 ### selected cell types to visualize
	colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
	NumCols = 4,                                 ### number of columns in the figure panel
        pointSize = 1.0)                             ### point size in ggplot2 scatterplot  
print(p2)


## VISUALIZE THE PROPORTION FOR TWO CELL TYPES

### visualize the spatial distribution of two cell types on the same plot
p3 <- CARD.visualize.prop.2CT(
    proportion = CARD_obj@Proportion_CARD,                             ### Cell type proportion estimated by CARD
    spatial_location = CARD_obj@spatial_location,                      ### spatial location information
    ct2.visualize = c("Cancer_clone_A","Cancer_clone_B"),              ### two cell types you want to visualize
    colors = list(c("lightblue","lightyellow","red"),c("lightblue","lightyellow","black")))       ### two color scales                             
print(p3)


## VISUALIZE THE CELL TYPE PROPORTION CORRELATION
p4 <- CARD.visualize.Cor(CARD_obj@Proportion_CARD,colors = NULL) # if not provide, we will use the default colors
print(p4)


# REFINED SPATIAL MAP

## IPUTATION ON THE NEWLY GRIDED SPATIAL LOCATIONS
CARD_obj = CARD.imputation(CARD_obj,NumGrids = 2000,ineibor = 10,exclude = NULL)


## Visualize the newly grided spatial locations to see if the shape is correctly detected. If not, the user can provide the row names of the excluded spatial location data into the CARD.imputation function
location_imputation = cbind.data.frame(x=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",1)),
	y=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",2)))
rownames(location_imputation) = rownames(CARD_obj@refined_prop)

library(ggplot2)

p5 <- ggplot(
    location_imputation, 
    aes(x = x, y = y)) + 
    geom_point(shape=22,color = "#7dc7f5")+
    theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.position="bottom",
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(colour = "grey89", fill=NA, size=0.5)
    )
print(p5)

## Visualize the cell type proportion at an enhanced resolution
p6 <- CARD.visualize.prop(
	proportion = CARD_obj@refined_prop,                         
	spatial_location = location_imputation,            
	ct.visualize = ct.visualize,                    
	colors = c("lightblue","lightyellow","red"),    
	NumCols = 4
    )                                  
print(p6)

## Visualize the marker gene expression at an enhanced resolution 
p7 <- CARD.visualize.gene(
	spatial_expression = CARD_obj@refined_expression,
	spatial_location = location_imputation,
	gene.visualize = c("Tm4sf1","S100a4","Tff3","Apol1","Crisp3","CD248"),
	colors = NULL,
	NumCols = 6
    )
print(p7)



# EXTENSION OF CARD IN A REFERENCE-FREE VERSION: CARDFREE

## CREATE AN CARDFREE OBJECT
### load the marker gene list
CARDfree_obj = createCARDfreeObject(
	markerList = markerList,
	spatial_count = spatial_count,
	spatial_location = spatial_location,
	minCountGene = 100,
	minCountSpot =5) 

## DECONVOLUTION USING CARDFREE
### deconvolution using CARDfree
CARDfree_obj = CARD_refFree(CARDfree_obj)


## VISUALIZATION OF THE RESULT OF CARDFREE
colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4","#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")
### In order to maximumply match with the original results of CARD, we order the colors to generally match with the results infered by CARD
CARDfree_obj@Proportion_CARD = CARDfree_obj@Proportion_CARD[,c(8,10,14,2,1,6,12,18,7,13,20,19,16,17,11,15,4,9,3,5)]
colnames(CARDfree_obj@Proportion_CARD) = paste0("CT",1:20)
p9 <- CARD.visualize.pie(
    CARDfree_obj@Proportion_CARD,CARDfree_obj@spatial_location,
    colors = colors
    )
print(p9)


# EXTENSION OF CARD FOR SINGLE CELL RESOLUTION MAPPING

####' Note that here the shapeSpot is the user defined variable 
####' which indicates the capturing area of single cells. 
####' Details see above.
scMapping = CARD_SCMapping(CARD_obj,shapeSpot="Square",numCell=20,ncore=10)
print(scMapping)


### spatial location info and expression count of the single cell resolution data
library(SingleCellExperiment)
MapCellCords = as.data.frame(colData(scMapping))
count_SC = assays(scMapping)$counts


library(ggplot2)
df = MapCellCords
colors = c("#8DD3C7","#CFECBB","#F4F4B9","#CFCCCF","#D1A7B9","#E9D3DE","#F4867C","#C0979F",
	"#D5CFD6","#86B1CD","#CEB28B","#EDBC63","#C59CC5","#C09CBF","#C2D567","#C9DAC3","#E1EBA0",
	"#FFED6F","#CDD796","#F8CDDE")
p10 = ggplot(df, aes(x = x, y = y, colour = CT)) + 
    geom_point(size = 3.0) +
    scale_colour_manual(values =  colors) +
    #facet_wrap(~Method,ncol = 2,nrow = 3) + 
        theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
            panel.background = element_rect(colour = "white", fill="white"),
            plot.background = element_rect(colour = "white", fill="white"),
    legend.position="bottom",
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
    axis.text =element_blank(),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 13,face="bold"),
    legend.text=element_text(size = 12),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(0.45, 'cm'),
    strip.text = element_text(size = 15,face="bold"))+
                                guides(color=guide_legend(title="Cell Type"))
print(p10)


#############

library("celldex")
ref <- ImmGenData(ensembl = TRUE, cell.ont = c("all", "nonna", "none"))


install_version(
  "BiocFileCache",
  version = 2.11.1,
  dependencies = NA,
  upgrade = c("default", "ask", "always", "never"),
  force = FALSE,
  quiet = FALSE,
  build = FALSE,
  build_opts = c("--no-resave-data", "--no-manual", "--no-build-vignettes"),
  build_manual = FALSE,
  build_vignettes = FALSE,
  repos = getOption("repos"),
  type = "source",
  ...
)

install_version("BiocFileCache", version = 2.10, repos = getOption("repos"), type = getOption("pkgType"), quiet = FALSE)


packageVersion("BiocFileCache")
packageVersion("dplyr")

remove.packages("BiocFileCache")

path_cellType <- "/mnt/d4t/SPATIAL/GeneSmart/CellType_annotated/"
ref_matrix <- readRDS(paste0(path_cellType,"GSE131907_Lung_Cancer_normalized_log2TPM_matrix.rds"))
View(ref_matrix)
