# Valiante.R
###########################################################################################################################- 

# load libraries ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(Seurat)
  library(SeuratData)
  library(Azimuth)
  library(patchwork)
  library(ggsci)
  library(ggpubr)
  library(gridExtra)
  library(svglite)
})
options(stringsAsFactors=FALSE)

#colors----
stallion = c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B","#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4","#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8","#6E4B9E","#0C727C","#7E1416","#D8A767","#3D3D3D")

# load data ----
matrix = ReadMtx(mtx = "./data/Valiante/matrix.mtx.gz", features = "./data/Valiante/features.tsv.gz", cells = "./data/Valiante/barcodes.tsv.gz")
# create seurat object 
seu_obj = CreateSeuratObject(matrix)
dim(seu_obj)
head(seu_obj)

# only consider cells with at least 200 detected genes and genes expressed in at least 5 cells
selected_c = WhichCells(seu_obj, expression = nFeature_RNA > 200)
selected_f = rownames(seu_obj)[rowSums(seu_obj) > 5]
seu_obj = subset(seu_obj, features = selected_f, cells = selected_c)
dim(seu_obj)
# basic qc check
seu_obj = PercentageFeatureSet(seu_obj, "^MT", col.name = "percent_mito")
seu_obj = PercentageFeatureSet(seu_obj, "^RP[SL]", col.name = "percent_ribo")
# violin plots 
feats = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo")
VlnPlot(seu_obj, features = feats, pt.size = 0, ncol = length(feats)) + NoLegend()

# run azimuth 
datasets = AvailableData()
seu_obj = RunAzimuth(seu_obj, reference = "humancortexref")
# plot clusters 
p1 = DimPlot(seu_obj, group.by = "predicted.subclass", pt.size = 1.1, label = TRUE, label.size = 4.5, repel = TRUE) + 
  scale_color_manual(values = stallion) + theme_classic2() + NoLegend() 

p4 = p2 + p3 + plot_layout(widths = c(1.75,1))

svglite(filename = "./figures/Valiante/grant_fig.svg", width = 5.5, height = 8)
p1 / p4 
dev.off()  

