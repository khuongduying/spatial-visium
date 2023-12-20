
# IMPORT LIBRARY
library(CARD)
library(MuSiC)
library(SingleCellExperiment)
library(SpatialExperiment)
library(readr)
library(Matrix)


# PATH
REF_CELLTYPE    <- "/mnt/rdisk/duydao/SPATIAL/GeneSmart/REF_CELLTYPE/"
rdata           <- "/mnt/rdisk/duydao/SPATIAL/GeneSmart/rdata/D1/"
fig_path        <- "/mnt/rdisk/duydao/SPATIAL/GeneSmart/figures/D1/"


# LOAD DATA

## SPATIAL DATA
SCE             <- readRDS(file=paste0(rdata," spatialD1_ready.rds")) 

## SPATIAL COUNT
spatial_count   <- counts(SCE) ## must be raw count

## SPATIAL LOCATION
spatial_location            <- colData(SCE)[c("row","col")]
### transformation
colnames(spatial_location)  <- c("x","y") ## transform to correct format
spatial_location            <- data.frame(spatial_location)
#spatial_location$x <- rev(spatial_location$x)


## SINGLE-CELL REF
## SINGLE-CELL REF COUNT

### load data
sc_count            <- readRDS(paste0(REF_CELLTYPE, "RNA_rawcounts_matrix.rds")) ## raw count

### transformation
### transform to a sparse matrix
sc_count            <- as(sc_count, "sparseMatrix")
sc_count            <- Matrix(data=sc_count, sparse=TRUE)

## SINGLE-CELL REF CELLTYPE ANNOTATED
### load reference database
sc_meta <- readRDS(paste0(REF_CELLTYPE, "refquery_final.rds"))

### View data
head(sc_meta[[]])

### Transformation
sc_meta[[]]["cellID"] <- rownames(sc_meta[[]])
sc_meta <- sc_meta[[]][c("cellID","predicted.celltypel2","Patient")]
colnames(sc_meta) <- c("cellID","cellType", "sampleInfo")

# CREATE AN CARD OBJECT
CARD_obj = createCARDObject(
	sc_count            = sc_count,
	sc_meta             = sc_meta,
	spatial_count       = spatial_count,
	spatial_location    = spatial_location,
	ct.varname          = "cellType",
	ct.select           = unique(sc_meta$cellType),
	sample.varname      = "sampleInfo",
	minCountGene        = 100,
	minCountSpot        = 5
    ) 


# DECONVOLUTION USING CARD

## DECONVOLUTION
CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
print(CARD_obj@Proportion_CARD[1:2,])


## VISUALIZE THE PROPORTION OF EACH CELL TYPE

### set the colors 
random_colors <- replicate(27, paste0("#", sprintf("%06x", sample(0:16777215, 1))))

### visualize
p1 <- CARD.visualize.pie(
	proportion          = CARD_obj@Proportion_CARD,
	spatial_location    = CARD_obj@spatial_location, 
    colors              = random_colors, 
    radius              = 0.4 ### You can choose radius = NULL or your own radius number
    ) + coord_flip() + scale_x_reverse()
print(p1)

### save figure
ggsave(
    paste0(fig_path,"CARD_proportion_D1.png"), 
    p1, 
    width   = 3000, 
    height  = 3000, 
    unit    = "px", 
    dpi     = 300
    )


### select the cell type that we are interested
ct.visualize = unique(sc_meta$cellType)
ct.visualize = ct.visualize[is.na(ct.visualize)==FALSE]
#ct.visualize = c("Exhausted CD8+ T","Tumor ECs","Malignant cells","mo-Mac")
ct.visualize1 = ct.visualize[1:9]
ct.visualize2 = ct.visualize[10:18]
ct.visualize3 = ct.visualize[19:27]
### visualize the spatial distribution of the cell type proportion
p2 <- CARD.visualize.prop(
                proportion         = CARD_obj@Proportion_CARD,        
                spatial_location   = CARD_obj@spatial_location, 
                ct.visualize       = ct.visualize3,                 ### selected cell types to visualize
                colors             = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
                NumCols            = 3,                                 ### number of columns in the figure panel
                pointSize          = 1.1
                ) +
                coord_flip() +
                scale_x_reverse()                             ### point size in ggplot2 scatterplot  
print(p2)

### save fig
ggsave(paste0(fig_path,"CARD_celltype27_D1_3.png"), p2, width=3000, height=3000, unit="px", dpi=300)


## VISUALIZE THE PROPORTION FOR TWO CELL TYPES

### visualize the spatial distribution of two cell types on the same plot
ct2.visualize <- c("CD8+ Effector memory T cells","CXCL1 Cancer cells")
p3 <- CARD.visualize.prop.2CT(
    proportion          = CARD_obj@Proportion_CARD,                             ### Cell type proportion estimated by CARD
    spatial_location    = CARD_obj@spatial_location,                      ### spatial location information
    ct2.visualize       = ct2.visualize,              ### two cell types you want to visualize
    colors              = list(
                            c("lightblue","lightyellow","red"),
                            c("lightblue","lightyellow","black"))
    ) + coord_flip() + scale_x_reverse()  + scale_size_manual(20)                   
print(p3)

### save fig
ggsave(paste0(fig_path,"CARD_2celltype27_Cancer_w_TCD8_D1.png"), p3, width=3000, height=3000, unit="px", dpi=300)


## VISUALIZE THE CELL TYPE PROPORTION CORRELATION
p4 <- CARD.visualize.Cor(
    CARD_obj@Proportion_CARD,
    colors = NULL
    ) # if not provide, we will use the default colors
print(p4)

### save fig
ggsave(paste0(fig_path,"CARD_correlation27_D1.png"), p4, width=3000, height=3000, unit="px", dpi=300)


# REFINED SPATIAL MAP

## IPUTATION ON THE NEWLY GRIDED SPATIAL LOCATIONS
CARD_obj = CARD.imputation(CARD_obj,NumGrids = 2000,ineibor = 10,exclude = NULL)


## Visualize the newly grided spatial locations to see if the shape is correctly detected. If not, the user can provide the row names of the excluded spatial location data into the CARD.imputation function
location_imputation = cbind.data.frame(
    x=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",1)),
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
ct.visualize = ct.visualize

p6 <- CARD.visualize.prop(
	proportion          = CARD_obj@refined_prop,                         
	spatial_location    = location_imputation,            
	ct.visualize        = ct.visualize1,                    
	colors              = NULL,    
	NumCols             = 3,
	pointSize           = 1.4
    ) +  
    theme_bw() + 
    coord_flip() + 
    scale_x_reverse()
                                    
print(p6)

ggsave(paste0(fig_path,"CARD_celltype27_enhanced_D1.png"), p6, width=3000, height=3000, unit="px", dpi=300)

## Visualize the marker gene expression at an enhanced resolution 
p7 <- CARD.visualize.gene(
	spatial_expression  = CARD_obj@refined_expression,
	spatial_location    = location_imputation,
	gene.visualize      = c("KRT7","IL12B", "KRT16", "KRT17","TP63","NAPSA"),
	colors              = NULL,
	NumCols             = 3
    ) + coord_flip() + scale_x_reverse()
print(p7)

ggsave(paste0(fig_path,"CARD_marker_enhanced_D1.png"), p7, width=4000, height=3000, unit="px", dpi=300)


# EXTENSION OF CARD FOR SINGLE CELL RESOLUTION MAPPING

####' Note that here the shapeSpot is the user defined variable 
####' which indicates the capturing area of single cells. 
####' Details see above.
scMapping = CARD_SCMapping(CARD_obj,shapeSpot="Square",numCell=7,ncore=40)
print(scMapping)

saveRDS(scMapping, file = paste0(rdata,"scMapping_27celltypes_D1.rds"), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

readRDS(paste0(rdata,"scMapping_27_celltypes_D1.rds"))
### spatial location info and expression count of the single cell resolution data
library(SingleCellExperiment)

MapCellCords    = as.data.frame(colData(scMapping))
count_SC        = assays(scMapping)$counts


library(ggplot2)
df      = MapCellCords
colors  = c("#ed3e24","#C59CC5","blue","#C2D567","#C9DAC3","#E1EBA0",
	"#FFED6F","#CDD796","#F8CDDE","#8DD3C7","#CFECBB","#F4F4B9","#CFCCCF","#D1A7B9","#E9D3DE","#F4867C","#C0979F",
	"#D5CFD6","#86B1CD","#CEB28B")

random_colors <- replicate(43, paste0("#", sprintf("%06x", sample(0:16777215, 1))))

p10 = ggplot(df, aes(x = x, y = y, colour = CT)) + 
        geom_point(size = 1) +
        scale_colour_manual(values =  random_colors) +
        #coord_flip() +
        #facet_wrap(~Method,ncol = 2,nrow = 3) + 
        theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
            panel.background = element_rect(colour = "white", fill="white"),
            plot.background = element_rect(colour = "white", fill="white"),
            legend.position="right",
            panel.border = element_rect(colour = "grey89", fill=NA, size=1),
            axis.text =element_blank(),
            axis.ticks =element_blank(),
            axis.title =element_blank(),
            legend.title=element_text(size = 13,face="bold"),
            legend.text=element_text(size = 13),
            legend.key = element_rect(colour = "transparent", fill = "white"),
            legend.key.size = unit(0.7, 'cm'),
            strip.text = element_text(size = 13,face="bold"))+
                                guides(color=guide_legend(title="Cell Type")) +
            coord_flip() + scale_x_reverse()
print(p10)


ggsave(paste0(fig_path,"CARD_celltype27_sc_res_D1_1.png"), p10, width=4000, height=3000, unit="px", dpi=300)


#############
