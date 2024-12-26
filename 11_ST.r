dev.off()
rm(list=ls())
dir()
gc()
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE) 
options(clusterProfiler.download.method = "wininet")
set.seed(1)
library(patchwork)
library(dplyr)
library(hdf5r)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(cowplot)
library(dplyr)
library(BayesSpace)
library(ArchR)
library(SeuratDisk)
library(patchwork)
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887") 
expr <- "/A4/raw_feature_bc_matrix/"
expr.mydata <- Seurat::Read10X(data.dir =  expr )
mydata <- Seurat::CreateSeuratObject(counts = expr.mydata, project = 'Posterior1', assay = 'Spatial')
mydata$slice <- "slice1"
mydata$region <- 'Posterior1' #命名
imgpath <- "/A4/spatial/"
img <- Seurat::Read10X_Image(image.dir = imgpath)
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = mydata)]
mydata[['image']] <- img
plot1 <- VlnPlot(mydata, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(mydata, features = "nCount_Spatial") + theme(legend.position = "right")
plot_grid(plot1, plot2)
plot1 <- VlnPlot(mydata, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(mydata, features = "nFeature_Spatial") + theme(legend.position = "right")
plot_grid(plot1, plot2)
plot1 <- VlnPlot(mydata, features = "nCount_Spatial", pt.size = 0.5) + NoLegend()
plot2 <- SpatialFeaturePlot(mydata, features = "nCount_Spatial") + theme(legend.position = "right")
plot_grid(plot1, plot2)
mydata <- SCTransform(mydata, assay = "Spatial", verbose = FALSE)
# rerun normalization to store sctransform residuals for all genes
mydata <- SCTransform(mydata, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
# also run standard log normalization for comparison
mydata <- NormalizeData(mydata, verbose = FALSE, assay = "Spatial") 
mydata <- GroupCorrelation(mydata, group.assay = "Spatial", assay = "Spatial", slot = "data", do.plot = FALSE)
mydata <- GroupCorrelation(mydata, group.assay = "Spatial", assay = "SCT", slot = "scale.data", do.plot = FALSE)
p1 <- GroupCorrelationPlot(mydata, assay = "Spatial", cor = "nCount_Spatial_cor") + ggtitle("Log Normalization") + 
  theme(plot.title = element_text(hjust = 0.5))
p2 <- GroupCorrelationPlot(mydata, assay = "SCT", cor = "nCount_Spatial_cor") + ggtitle("SCTransform Normalization") + 
  theme(plot.title = element_text(hjust = 0.5))
p3 <- plot_grid(p1, p2)
p3
SpatialFeaturePlot(mydata, features = c("THBS2", "COL1A1","EPCAM"))
















