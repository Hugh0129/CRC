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

colon_stromal = readRDS("./colon_stromal.rds")

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
	# print(pdf(paste0("./pca", sample_name, ".pdf")))
	# DimPlot(colon, reduction = "pca")
	# dev.off()
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

# Idents(colon) = "seurat_clusters"
DotPlot(colon, 
        features =c("PECAM1", "VWF",
                    "BTNL9",'LNX1',"ATP2A3",# arterial
                    "SERPINE1",'COL4A1','COL4A2','COL15A1',#capillary
                    "LYVE1","PROX1","RELN","EFNA5", 'MMRN1',# lymphatic
                    "CPE","ADGRG6","STXBP6",'IL1R1', # venous
                    "MKI67",'CEP55', "TOP2A",'NUSAP1'# 'UBE2C',"TYMS",  #'POF1B','RNF43','RBFOX1' #proliferating 
                   ), 
        # cols = c("blue", "red"), 
        dot.scale = 8, 
        # split.by = "tissue"
        ) + RotatedAxis()

colon@meta.data$celltype = "NA"
lymphatic <- c(4,7,9)
proliferating <- c(15)
venous <- c(2,10,13)
capillary <- c(0,1,3,6,12,14)
arterial <- c(5,8,11)
colon@meta.data[which(colon@meta.data$seurat_clusters %in% lymphatic),'celltype'] <- "lymphatic"
colon@meta.data[which(colon@meta.data$seurat_clusters %in% proliferating),'celltype'] <- "proliferating"
colon@meta.data[which(colon@meta.data$seurat_clusters %in% venous),'celltype'] <- "venous" # 静脉
colon@meta.data[which(colon@meta.data$seurat_clusters %in% capillary),'celltype'] <- "capillary" # 毛细血管
colon@meta.data[which(colon@meta.data$seurat_clusters %in% arterial),'celltype'] <- "arterial" # 动脉
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
           # split.by = "DiseaseState",
            pt.size = 1,
            label = F)

# Idents(colon) = "seurat_clusters"
dotplot_endo = DotPlot(colon_endo, 
        features =c('EFNB2',"ATP2A3",# arterial
                    "ADGRG6", 'IL1R1', # venous
                    "RELN","EFNA5", # lymphatic
                    "SERPINE1",'COL4A2',#capillary'COL4A1',
                    "MKI67",'NUSAP1'# 'UBE2C',"TYMS",  #'POF1B','RNF43','RBFOX1' #proliferating  ???
                    #'PTPRC'#myogenic
                    #'CSPG4','CD34','MCAM' # pericyte
                    # 'PECAM1','CD34' # mature endo
                    
                   ), 
        # cols = c("blue", "red"), 
        dot.scale = 8, 
        # split.by = "tissue"
        ) + RotatedAxis()
dotplot_endo

cell.prop<-as.data.frame(prop.table(table(colon_endo@meta.data$celltype,colon_endo$DiseaseState)))

colnames(cell.prop)<-c("cluster","sample","proportion")

library("ggplot2")

p = ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  labs(x = '', y = 'Relative Abundance(%)') +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13)) +
  theme(legend.text = element_text(size = 11))

# colors
ccc = c('#e6194b', '#3cb44b',  '#4363d8', '#f58231', '#911eb4', 
        '#46f0f0', '#f032e6', '#bcf60c',  '#008080', 
        '#800000',  '#808000',  '#000075', '#808080')
cell_type_cols <- c(brewer.pal(9, "Set1"), 
                    # brewer.pal(9, "Set1"), 
                    "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F",
                    "#FF6A6A","#7FFFD4", "#AB82FF","#90EE90","#00CD00","#008B8B",
                    "#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00")  

#cell_type_cols <-paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE)
#cell_type_cols <-paletteDiscrete(values = unique(colon@meta.data$celltype1), set = "stallion", reverse = FALSE)
cell_type_cols <-paletteDiscrete(values = unique(colon_endo@meta.data$celltype), set = "stallion", reverse = FALSE)

p <- p + scale_fill_manual(values = (cell_type_cols)) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank())+
  guides(color = guide_legend(ncol = 1, byrow = TRUE,reverse = T))+
  theme(axis.title.y = element_text(face = 'plain',color = 'black',size = 16),
        axis.title.x = element_text(face = 'plain',color = 'black',size = 10),
        axis.text.y = element_text(face = 'plain',color = 'black',size = 16),
        axis.text.x = element_text(face = 'plain',color = 'black',
                                   size = 16,angle = 90,vjust = 0.5, hjust=0), #x轴标签偏转45°，并下降0.5
        axis.ticks.length=unit(.1,"lines"),#设置刻度线的高度
        axis.ticks.x = element_line(size=0.05, colour = "black"),
        #axis.ticks.margin=unit(.4,"cm"),#设置刻度数字与刻度线的距离
        panel.border = element_blank(),
        axis.line = element_line(size=0.1, colour = "black"),
        #axis.ticks.x.bottom =  = element_line(size = 0.5),
        panel.grid = element_blank(),
        #scale_fill_distiller(palette = "Spectral"),
        #scale_fill_brewer(palette = 'Paired'),
        #scale_fill_manual(values=ccc),
        legend.position = 'right',
        legend.key.width = unit(0.2,'cm'),
        legend.key.height = unit(0.2,'cm'),#定义图例中色块的高度
        legend.text = element_text(color = 'black',size = 16))
p








