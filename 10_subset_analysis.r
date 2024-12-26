# Endothelial
getwd()
rm(list=ls())
dir()
gc()
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ArchR)
library(viridis)
library(DoubletFinder)
library(Rcpp)
library(harmony)
library(future)
library(Matrix)
library(Rmagic)
library(ggpubr)
set.seed(1)



colon = readRDS("./colon_Endothelial_all_samples_initial.rds")

# Define variables
sample_name <- "all_samples" # used as label for saving plots
execute_steps <- c(1,2,3)
vln_plot <- function(features, save_name){
	pdf(save_name, width = 20, onefile=F)
	print(VlnPlot(colon, features = features, group.by = "orig.ident", pt.size = 0, cols = c(rep("#D51F26",100)))+
	geom_boxplot(outlier.shape = NA, alpha = 0.6)+theme_ArchR()+theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))+
	scale_x_discrete(labels=paste0(data.frame(table(colon@meta.data$orig.ident))$Var1, "\n n = ",  data.frame(table(colon@meta.data$orig.ident))$Freq)))
	dev.off()
}

normalize_and_dim_reduce <- function (colon, sample_name){
	# stanfard seurat pipeline
	colon <- NormalizeData(colon, normalization.method = "LogNormalize", scale.factor = 10000)
	colon <- FindVariableFeatures(colon, selection.method = "vst", nfeatures = 2000)
	all.genes <- rownames(colon)
	colon <- ScaleData(colon, features = all.genes)
	colon <- RunPCA(colon, features = VariableFeatures(object = colon))
	print(pdf(paste0("./pca", sample_name, ".pdf")))
	DimPlot(colon, reduction = "pca")
	dev.off()
	colon <- FindNeighbors(colon, dims = 1:30)
	colon <- FindClusters(colon, resolution = 1.0)
	colon <- RunUMAP(colon, dims = 1:30)
	return(colon)
}
plotUMAPandRunHarmony <- function(colon, run_harmony, version){
	 plot umaps
	 pdf(paste0("./UMAPclustering", version, ".pdf"))
	 print(DimPlot(colon, reduction = "umap", cols = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE)))
	 dev.off()
	 pdf(paste0("./UMAP_samples", version, ".pdf"), width = 12)
	 print(DimPlot(colon, reduction = "umap", group.by = "orig.ident", cols = paletteDiscrete(values = unique(colon@meta.data$orig.ident), set = "stallion", reverse = FALSE)))
	 dev.off()
	 pdf(paste0("./UMAP_disease_state", version, ".pdf"), width = 6.5, onefile=F)
	 print(DimPlot(colon, reduction = "umap", group.by = "DiseaseState",
	 	cols = c("#D51F26", "#272E6A", "#208A42", "#89288F" ,"#F47D2B" ,"#FEE500" ,"#8A9FD1", "#C06CAB","#90D5E4", "#89C75F")) + theme_ArchR())
	 dev.off()
	
	# run harmony
	if (run_harmony){
		library(harmony)
		colon <- RunHarmony(colon, "orig.ident") # will use PCA
		colon <- RunUMAP(colon, dims = 1:30, reduction = "harmony", reduction.name = "umapharmony")
		 pdf(paste0("./UMAPharmony_samples.pdf"), width = 12)
		 print(DimPlot(colon, reduction = "umapharmony", group.by = "orig.ident", cols = paletteDiscrete(values = unique(colon@meta.data$orig.ident), set = "stallion", reverse = FALSE)))#, cols = (ArchRPalettes$stallion))
		 dev.off()

		 pdf(paste0("./UMAPharmony_diseasestate.pdf"), width = 12)
		 print(DimPlot(colon, reduction = "umapharmony", group.by = "DiseaseState", cols = paletteDiscrete(values = unique(colon@meta.data$DiseaseState), set = "stallion", reverse = FALSE)))#, cols = (ArchRPalettes$stallion))
		 dev.off()

		 pdf(paste0("./UMAPharmony_clustering.pdf"), width = 12)
		 print(DimPlot(colon, reduction = "umapharmony", cols = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE)))
		 dev.off()

		colon <- FindNeighbors(colon, reduction = "harmony", dims = 1:30)
		colon <- FindClusters(colon, resolution = 1.0)

		 pdf(paste0("./UMAPharmony_clustering_new.pdf"), width = 12)
		 print(DimPlot(colon, reduction = "umapharmony", cols = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE)))
		 dev.off()
	}
	return(colon)
}

# 1) Normalize and scale data
if (1 %in% execute_steps){
	colon <- normalize_and_dim_reduce(colon, sample_name)
}
# 2) Plot UMAPs and run harmony
if (2 %in% execute_steps){
	colon <- plotUMAPandRunHarmony(colon, TRUE, version = "initial")
}
# 3) Remove bad clusters and redo analysis
if (3 %in% execute_steps){
	# Save ids of cells to remove
	bad_clusters <- bad_clusters <- c(2,5,11,12)
	# Subset
	colon_new <- DietSeurat(subset(colon, subset = seurat_clusters %ni% bad_clusters))
	# Redo the initial analysis steps
	colon_new <- normalize_and_dim_reduce(colon_new, sample_name)
	colon_new <- FindClusters(colon_new, resolution = 2.0)
	# colon <- colon_new
}
DimPlot(colon, reduction = "umap",label = T,
        cols = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE))

colon@meta.data$celltype = "NA"
lymphatic <- c(4,7,9)
proliferating <- c(15)
venous <- c(2,10,13)
capillary <- c(0,1,3,6,12,14)
arterial <- c(5,8,11)
colon@meta.data[which(colon@meta.data$seurat_clusters %in% lymphatic),'celltype'] <- "lymphatic"
colon@meta.data[which(colon@meta.data$seurat_clusters %in% proliferating),'celltype'] <- "proliferating"
colon@meta.data[which(colon@meta.data$seurat_clusters %in% venous),'celltype'] <- "venous" 
colon@meta.data[which(colon@meta.data$seurat_clusters %in% capillary),'celltype'] <- "capillary" 
colon@meta.data[which(colon@meta.data$seurat_clusters %in% arterial),'celltype'] <- "arterial"
Idents(colon) = "celltype"
DimPlot(colon, reduction = "umap",label = F,
        cols = paletteDiscrete(values = unique(colon@meta.data$celltype), set = "stallion", reverse = FALSE))
DimPlot(colon, reduction = "umap",label = F,
        cols = paletteDiscrete(values = unique(colon@meta.data$celltype), set = "stallion", reverse = FALSE),
       split.by = "DiseaseState")

FeaturePlot(colon,
            reduction = "umap",
            features = c('CD36'),
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            pt.size = 1,
            label = F)

dotplot_endo = DotPlot(colon, 
        features =c('EFNB2',"ATP2A3",# arterial
                    "ADGRG6", 'IL1R1', # venous
                    "RELN","EFNA5", # lymphatic
                    "SERPINE1",'COL4A2',#capillary
                    "MKI67",'NUSAP1'#proliferating 
                   ), 
        dot.scale = 8, 
        ) + RotatedAxis()
dotplot_endo

cell.prop<-as.data.frame(prop.table(table(colon@meta.data$celltype,colon$DiseaseState)))
colnames(cell.prop)<-c("cluster","sample","proportion")
p = ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  labs(x = '', y = 'Relative Abundance(%)') +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13)) +
  theme(legend.text = element_text(size = 11))


cell_type_cols <-paletteDiscrete(values = unique(colon@meta.data$celltype), set = "stallion", reverse = FALSE)

p <- p + scale_fill_manual(values = (cell_type_cols)) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank())+
  guides(color = guide_legend(ncol = 1, byrow = TRUE,reverse = T))+
  theme(axis.title.y = element_text(face = 'plain',color = 'black',size = 16),
        axis.title.x = element_text(face = 'plain',color = 'black',size = 10),
        axis.text.y = element_text(face = 'plain',color = 'black',size = 16),
        axis.text.x = element_text(face = 'plain',color = 'black',
                                   size = 16,angle = 90,vjust = 0.5, hjust=0),
        axis.ticks.length=unit(.1,"lines"),
        axis.ticks.x = element_line(size=0.05, colour = "black"),
        panel.border = element_blank(),
        axis.line = element_line(size=0.1, colour = "black"),
        panel.grid = element_blank(),
        legend.position = 'right',
        legend.key.width = unit(0.2,'cm'),
        legend.key.height = unit(0.2,'cm'),
        legend.text = element_text(color = 'black',size = 16))
p

#=============================================================================================================================
# Fibroblasts
getwd()
rm(list=ls())
dir()
gc()
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ArchR)
library(viridis)
library(DoubletFinder)
library(Rcpp)
library(harmony)
library(future)
library(Matrix)
library(Rmagic)
library(ggpubr)
set.seed(1)
colon = readRDS("./colon_Fibroblasts_all_samples_initial.rds")
# Define variables
sample_name <- "all_samples" # used as label for saving plots
execute_steps <- c(1,2,3,4)
vln_plot <- function(features, save_name){
	pdf(save_name, width = 20, onefile=F)
	print(VlnPlot(colon, features = features, group.by = "orig.ident", pt.size = 0, cols = c(rep("#D51F26",100)))+
	geom_boxplot(outlier.shape = NA, alpha = 0.6)+theme_ArchR()+theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))+
	scale_x_discrete(labels=paste0(data.frame(table(colon@meta.data$orig.ident))$Var1, "\n n = ",  data.frame(table(colon@meta.data$orig.ident))$Freq)))
	dev.off()
}
normalize_and_dim_reduce <- function (colon, sample_name){
	# stanfard seurat pipeline
	colon <- NormalizeData(colon, normalization.method = "LogNormalize", scale.factor = 10000)
	colon <- FindVariableFeatures(colon, selection.method = "vst", nfeatures = 2000)
	all.genes <- rownames(colon)
	colon <- ScaleData(colon, features = all.genes)
	colon <- RunPCA(colon, features = VariableFeatures(object = colon))
	print(pdf(paste0("./pca", sample_name, ".pdf")))
	DimPlot(colon, reduction = "pca")
	dev.off()
	colon <- FindNeighbors(colon, dims = 1:30)
	colon <- FindClusters(colon, resolution = 2.0)
	colon <- RunUMAP(colon, dims = 1:30)
	return(colon)
}

plotUMAPandRunHarmony <- function(colon, run_harmony, version){
	# plot umaps
	pdf(paste0("./UMAPclustering", version, ".pdf"))
	print(DimPlot(colon, reduction = "umap", cols = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE)))
	dev.off()
	pdf(paste0("./UMAP_samples", version, ".pdf"), width = 12)
	print(DimPlot(colon, reduction = "umap", group.by = "orig.ident", cols = paletteDiscrete(values = unique(colon@meta.data$orig.ident), set = "stallion", reverse = FALSE)))
	dev.off()
	pdf(paste0("./UMAP_disease_state", version, ".pdf"), width = 6.5, onefile=F)
	print(DimPlot(colon, reduction = "umap", group.by = "DiseaseState",
		cols = c("#D51F26", "#272E6A", "#208A42", "#89288F" ,"#F47D2B" ,"#FEE500" ,"#8A9FD1", "#C06CAB","#90D5E4", "#89C75F")) + theme_ArchR())
	dev.off()
	# run harmony
	if (run_harmony){
		library(harmony)
		colon <- RunHarmony(colon, "orig.ident") # will use PCA
		colon <- RunUMAP(colon, dims = 1:30, reduction = "harmony", reduction.name = "umapharmony")
		pdf(paste0("./UMAPharmony_samples.pdf"), width = 12)
		print(DimPlot(colon, reduction = "umapharmony", group.by = "orig.ident", cols = paletteDiscrete(values = unique(colon@meta.data$orig.ident), set = "stallion", reverse = FALSE)))
		dev.off()
		pdf(paste0("./UMAPharmony_diseasestate.pdf"), width = 12)
		print(DimPlot(colon, reduction = "umapharmony", group.by = "DiseaseState", cols = paletteDiscrete(values = unique(colon@meta.data$DiseaseState), set = "stallion", reverse = FALSE)))
		dev.off()
		pdf(paste0("./UMAPharmony_clustering.pdf"), width = 12)
		print(DimPlot(colon, reduction = "umapharmony", cols = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE)))
		dev.off()
		colon <- FindNeighbors(colon, reduction = "harmony", dims = 1:30)
		colon <- FindClusters(colon, resolution = 2.0)
		pdf(paste0("./UMAPharmony_clustering_new.pdf"), width = 12)
		print(DimPlot(colon, reduction = "umapharmony", cols = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE)))
		dev.off()
	}
	return(colon)
}
# 1) Normalize and scale data
if (1 %in% execute_steps){
	colon <- normalize_and_dim_reduce(colon, sample_name)
}
# 2) Plot UMAPs and run harmony
if (2 %in% execute_steps){
	colon <- plotUMAPandRunHarmony(colon, TRUE, version = "initial")
}
# 3) Remove bad clusters and redo analysis
if (3 %in% execute_steps){
	# Save ids of cells to remove
	bad_clusters <- bad_clusters <- c(8,10,13,16,24,25,26,27)
	# Subset
	colon_new <- DietSeurat(subset(colon, subset = seurat_clusters %ni% bad_clusters))
	# Redo the initial analysis steps
	colon_new <- normalize_and_dim_reduce(colon_new, sample_name)
	colon_new <- FindClusters(colon_new, resolution = 2.0)
	colon <- colon_new
}
# 4) Remove bad clusters and redo analysis
if (4 %in% execute_steps){
	# Save ids of cells to remove
	bad_clusters <- bad_clusters <- c(9,14)
	# Subset
	colon_new <- DietSeurat(subset(colon, subset = seurat_clusters %ni% bad_clusters))
	# Redo the initial analysis steps
	colon_new <- normalize_and_dim_reduce(colon_new, sample_name)
	colon_new <- FindClusters(colon_new, resolution = 2.0)
	colon <- colon_new
}
DimPlot(colon, reduction = "umap",label = T,
        cols = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE))
colon@meta.data$celltype = "NA"
NF <- c(0,1,3,6,8,10,12,13,14,16,18,19,22,26)
CAF_myo <- c(7,11,20,27)
CAF_infla <- c(2,4,5,9,15,21,23,25)
CAF_adi <- c(24)
CAF_PN <- c(17,28)
colon@meta.data[which(colon@meta.data$seurat_clusters %in% NF),'celltype'] <- "NF"
colon@meta.data[which(colon@meta.data$seurat_clusters %in% CAF_myo),'celltype'] <- "CAF_myo"
colon@meta.data[which(colon@meta.data$seurat_clusters %in% CAF_infla),'celltype'] <- "CAF_infla"
colon@meta.data[which(colon@meta.data$seurat_clusters %in% CAF_adi),'celltype'] <- "CAF_adi"
colon@meta.data[which(colon@meta.data$seurat_clusters %in% CAF_PN),'celltype'] <- "CAF_PN"
Idents(colon) = "celltype"
DimPlot(colon, reduction = "umap",label = F,
        cols = paletteDiscrete(values = unique(colon@meta.data$celltype), set = "stallion", reverse = FALSE))
TH = FeaturePlot(colon,
            reduction = "umap",
            features = c("THBS2"),
            sort.cell = TRUE,
            pt.size = 1,
            raster=FALSE,
            label = F)
TH
CAFs_clusters = c("CAF_adi","CAF_infla","CAF_myo","CAF_PN")
colon_other = subset(colon, subset = celltype %in% CAFs_clusters)
dotplot_caf = DotPlot(colon_other, 
        features =c("COL1A1","COL1A2","COL3A1", "MYH11","TPM1",
                    "GSN",'NRXN1',"GPM6B"
                   ), 
        # cols = c("blue", "red"), 
        dot.scale = 8, 
        # split.by = "tissue"
        ) + RotatedAxis()
dotplot_caf
cell.prop<-as.data.frame(prop.table(table(colon@meta.data$celltype,colon$DiseaseState)))
colnames(cell.prop)<-c("cluster","sample","proportion")
p = ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  labs(x = '', y = 'Relative Abundance(%)') +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13)) +
  theme(legend.text = element_text(size = 11))
cell_type_cols <-paletteDiscrete(values = unique(colon@meta.data$celltype), set = "stallion", reverse = FALSE)
p <- p + scale_fill_manual(values = (cell_type_cols)) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank())+
  guides(color = guide_legend(ncol = 1, byrow = TRUE,reverse = T))+
  theme(axis.title.y = element_text(face = 'plain',color = 'black',size = 16),
        axis.title.x = element_text(face = 'plain',color = 'black',size = 10),
        axis.text.y = element_text(face = 'plain',color = 'black',size = 16),
        axis.text.x = element_text(face = 'plain',color = 'black',
                                   size = 16,angle = 90,vjust = 0.5, hjust=0), 
        axis.ticks.length=unit(.1,"lines"),
        axis.ticks.x = element_line(size=0.05, colour = "black"),
        panel.border = element_blank(),
        axis.line = element_line(size=0.1, colour = "black"),
        panel.grid = element_blank(),
        legend.position = 'right',
        legend.key.width = unit(0.2,'cm'),
        legend.key.height = unit(0.2,'cm'),
        legend.text = element_text(color = 'black',size = 16))
p


































