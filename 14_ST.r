# normal colon 
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
expr <- "./A4/raw_feature_bc_matrix/"
expr.mydata <- Seurat::Read10X(data.dir =  expr )
mydata <- Seurat::CreateSeuratObject(counts = expr.mydata, project = 'Posterior1', assay = 'Spatial')
mydata$slice <- "slice1"
mydata$region <- 'Posterior1' #命名
imgpath <- "./A4/spatial/"
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

#==================================================================================================================================
# CRC
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
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887") 
CRC<-Load10X_Spatial(data.dir = 'data', assay = "Spatial",slice = "slice1", filter.matrix = TRUE, to.upper = TRUE )
CRC@images[["slice1"]]@coordinates[["tissue"]] <- as.integer(CRC@images[["slice1"]]@coordinates[["tissue"]])
CRC@images[["slice1"]]@coordinates[["row"]] <- as.integer(CRC@images[["slice1"]]@coordinates[["row"]])
CRC@images[["slice1"]]@coordinates[["col"]] <- as.integer(CRC@images[["slice1"]]@coordinates[["col"]])
CRC@images[["slice1"]]@coordinates[["imagerow"]] <- as.integer(CRC@images[["slice1"]]@coordinates[["imagerow"]])
CRC@images[["slice1"]]@coordinates[["imagecol"]] <- as.integer(CRC@images[["slice1"]]@coordinates[["imagecol"]])
plot1 <- VlnPlot(CRC, features = "nCount_Spatial", pt.size = 0.5) + NoLegend()
plot2 <- SpatialFeaturePlot(CRC, features = "nCount_Spatial") + theme(legend.position = "right")
plot_grid(plot1, plot2)
CRC <- SCTransform(CRC, assay = "Spatial", verbose = FALSE)
# rerun normalization to store sctransform residuals for all genes
CRC <- SCTransform(CRC, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
# also run standard log normalization for comparison
CRC <- NormalizeData(CRC, verbose = FALSE, assay = "Spatial")
# Computes the correlation of the log normalized data and sctransform residuals with the number of UMIs
CRC <- GroupCorrelation(CRC, group.assay = "Spatial", assay = "Spatial", slot = "data", do.plot = FALSE)
CRC <- GroupCorrelation(CRC, group.assay = "Spatial", assay = "SCT", slot = "scale.data", do.plot = FALSE)
p1 <- GroupCorrelationPlot(CRC, assay = "Spatial", cor = "nCount_Spatial_cor") + ggtitle("Log Normalization") + 
  theme(plot.title = element_text(hjust = 0.5))
p2 <- GroupCorrelationPlot(CRC, assay = "SCT", cor = "nCount_Spatial_cor") + ggtitle("SCTransform Normalization") + 
  theme(plot.title = element_text(hjust = 0.5))
p3 <- plot_grid(p1, p2)
p3
CRC <- RunPCA(CRC, assay = "SCT", verbose = FALSE)
CRC <- FindNeighbors(CRC, reduction = "pca", dims = 1:30)
CRC <- FindClusters(CRC, verbose = FALSE)
CRC <- RunUMAP(CRC, reduction = "pca", dims = 1:30)
p1 <- DimPlot(CRC, reduction = "umap", label = TRUE,
              cols = paletteDiscrete(values = unique(CRC@meta.data$seurat_clusters)))
p2 <- SpatialDimPlot(CRC, label = TRUE, label.size = 3,
                     cols = paletteDiscrete(values = unique(CRC@meta.data$seurat_clusters)))
plot_grid(p1, p2)
SpatialFeaturePlot(CRC, features = c("THBS2", "EPCAM"))
FeaturePlot(CRC,
            reduction = "umap",
            features = c("THBS2"),
            sort.cell = TRUE,
            pt.size = 1,
            label = TRUE)
CRC@meta.data$Celltype = "NA"
CRC@meta.data$Celltype[which(CRC@meta.data$seurat_clusters%in%c(9,10,11,14))] <- 'THBS2_CAF'
CRC@meta.data$Celltype[which(CRC@meta.data$seurat_clusters%in%c(1))] <- 'CAF/Epi '
CRC@meta.data$Celltype[which(CRC@meta.data$seurat_clusters%in%c(16))] <- 'Myofib'
CRC@meta.data$Celltype[which(CRC@meta.data$seurat_clusters%in%c(0,2,3,4,5,8,13))] <- 'Epi'
CRC@meta.data$Celltype[which(CRC@meta.data$seurat_clusters%in%c(6))] <- 'SPP1_mac'
CRC@meta.data$Celltype[which(CRC@meta.data$seurat_clusters%in%c(17))] <- 'SPP1_mac/T/Epi'
CRC@meta.data$Celltype[which(CRC@meta.data$seurat_clusters%in%c(12))] <- 'Mac'
CRC@meta.data$Celltype[which(CRC@meta.data$seurat_clusters%in%c(7))] <- 'THBS2_CAF/T/B'
CRC@meta.data$Celltype[which(CRC@meta.data$seurat_clusters%in%c(15))] <- 'CAF/Endo'
Idents(CRC) = "Celltype"
p1 <- DimPlot(CRC, reduction = "umap", label = TRUE,
              cols = paletteDiscrete(values = unique(CRC@meta.data$Celltype)))
p2 <- SpatialDimPlot(CRC, label = TRUE, label.size = 3,
                     cols = paletteDiscrete(values = unique(CRC@meta.data$Celltype)))
plot_grid(p1, p2)
markers.to.plot <- c( "EPCAM","KRT8","CDH1", # epi
                      "PTPRC", # immune
                      "CD3D","CD3E", # T
                      "PECAM1","VWF", # endo
                      "COL1A1","COL3A1","MYLK","DES", "MYL9","THBS2", # fib
                      "CD14","CD68","SPP1", # mac/mon
                      "CPA3", # mast
                      "CD79A" # B
)
DotPlot(CRC, 
        features =markers.to.plot, 
        dot.scale = 8) + RotatedAxis()





