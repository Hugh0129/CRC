dev.off()
rm(list=ls())
dir()
gc()
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ArchR)
library(viridis)
library(DoubletFinder)
set.seed(1)
execute_steps <- c(1,2,3,4,5)
# Define variables
sample_name <- "all_samples"
individual_qc_and_dublet_plot_location <- "./initial_clustering/doublet_analysis/"
analysis_parent_folder <- "./initial_clustering/"
setwd(analysis_parent_folder)
path_to_metadata <- "hubmap_htan_metadata_atac_and_rna_final.csv"
# Define sets and locations of files for initial processing
scRNA_data_path <- "/data/"
scRNA_set_names <- c("all_samples_s12a", "all_samples_s12b", "all_samples_s5", "all_samples_s6", "all_samples_s7", "all_samples_s8", "all_samples_s3", "all_samples_s4", "all_samples_s9", "all_samples_s10", "all_samples_s11")
scRNA_sets <- 
  list(c("EP007","A014-C-201","A002-C-010","A001-C-207","A001-C-124"), 
			c("B004-A-204", "B004-A-104", "A002-C-106", "A001-C-223", "A002-C-204", "A014-C-111"),
			c("A015-C-203","A015-C-204","A014-C-040","A002-C-201","A002-C-203"), 
			c("B001-A-301","B001-A-401","A015-C-008","A015-C-208","A014-C-114","A014-C-043","A001-C-202", "A001-C-119"),
			c("A001-C-108","A002-C-021","A002-C-212","A002-C-205","A014-C-101","A014-C-108", "A001-C-104"),
			c("A015-C-005","A015-C-006","A015-C-106","A002-C-114","A015-C-104","A015-C-202"), 
			c("A001-C-014","A001-C-023","A002-C-016","A002-C-024","A014-C-001","A014-C-054","A015-C-002","A015-C-010"),
			c("A015-C-109","B004-A-004","A002-C-121-R0", "A002-C-010-R0"), 
			c("A001-C-007", "A001-C-203", "A002-C-116", "A002-C-121", "A008-E-008", "A008-E-015", "A010-E-018", "A010-E-023", "A014-C-008", "A014-C-052"),
			c("A015-C-001", "A018-E-013", "A018-E-020", "B001-A-406", "B001-A-501", "B004-A-008", "EP034", "EP072B", "EP091", "A022-E-022"),
			c("CRC1_8810", "CRC2_15564", "CRC3_11773")) 
# Define functions
seurat_standard_normalize_and_scale <- function(colon, cluster, cluster_resolution){
	# colon is seurat object, 
	colon <- NormalizeData(colon, normalization.method = "LogNormalize", scale.factor = 10000)
	colon <- FindVariableFeatures(colon, selection.method = "vst", nfeatures = 2000)
	all.genes <- rownames(colon)
	colon <- ScaleData(colon, features = all.genes)
	colon <- RunPCA(colon, features = VariableFeatures(object = colon))
	if (cluster){
		colon <- FindNeighbors(colon, dims = 1:20)
		colon <- FindClusters(colon, resolution = cluster_resolution)
	}
	colon <- RunUMAP(colon, dims = 1:20)
	return(colon)
}
make_seurat_object_and_doublet_removal <- function(data_directory, project_name){
	# function for basic seurat based qc and doubletfinder based doublet removal
	colon.data <- Read10X(data.dir = data_directory)
	currentSample <- CreateSeuratObject(counts = colon.data, project = project_name, min.cells = 3, min.features = 40)
	currentSample[["percent.mt"]] <- PercentageFeatureSet(currentSample, pattern = "^MT-")
	# qc plot-pre filtering
	pdf(paste0("./qc_plots_", project_name, "_prefiltered.pdf"))
	print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.05))
	dev.off()
	pdf(paste0("./qc_plots_", project_name, "_prefiltered_no_points.pdf"))
	print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0))
	dev.off()
	# filter everything to 400 unique genes/cell
	currentSample <- subset(currentSample, subset = nFeature_RNA > 400 & nFeature_RNA < 4000)
	# Normalize and make UMAP
	currentSample <- seurat_standard_normalize_and_scale(currentSample, FALSE)
	# Run doublet finder
	nExp_poi <- round(0.08*length(currentSample@meta.data$orig.ident)*length(currentSample@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
	seu_colon <- doubletFinder_v3(currentSample, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
	print(head(seu_colon@meta.data))
	# rename columns
	seu_colon$doublet.class <- seu_colon[[paste0("DF.classifications_0.25_0.09_",nExp_poi)]]
	seu_colon[[paste0("DF.classifications_0.25_0.09_",nExp_poi)]] <- NULL
	pann <- grep(pattern="^pANN", x=names(seu_colon@meta.data), value=TRUE)
	seu_colon$pANN <- seu_colon[[pann]]
	seu_colon[[pann]] <- NULL
	# plot pre and post doublet finder results
	pdf(paste0("./UMAP_pre_double_removal", project_name, ".pdf"))
	print(DimPlot(seu_colon, reduction = "umap", group.by = "doublet.class", cols = c("#D51F26", "#272E6A")))
	dev.off()
	seu_colon <- subset(seu_colon, subset = doublet.class != "Doublet")
	pdf(paste0("./UMAP_post_double_removal", project_name, ".pdf"))
	print(DimPlot(seu_colon, reduction = "umap", cols = c("#D51F26")))
	dev.off()
	# Remove extra stuff and return filtered Seurat object
	seu_colon <- DietSeurat(seu_colon, counts=TRUE, data=TRUE, scale.data=FALSE, assays="RNA")
	return(seu_colon)
}
seurat_qc_plots <- function(colon, sample_name){
	# Make some basic qc plots
	pdf(paste0("./seurat_nFeature_plots_", sample_name, ".pdf"), width = 40, height = 15)
	print(VlnPlot(colon, features = c("nFeature_RNA"), ncol = 1, pt.size = 0.2))
	dev.off()
	pdf(paste0("./seurat_nCount_plots_", sample_name, ".pdf"), width = 40, height = 15)
	print(VlnPlot(colon, features = c("nCount_RNA"), ncol = 1, pt.size = 0.2))
	dev.off()
	pdf(paste0("./seurat_pMT_plots_", sample_name, ".pdf"), width = 40, height = 15)
	print(VlnPlot(colon, features = c("percent.mt"), ncol = 1, pt.size = 0.2))
	dev.off()
}
# 1) Create seurat objects for individual 10x runs, run doblet finder and filter most likely doublets, merge into a seurat object containing all samples
if (1 %in% execute_steps){
	setwd(individual_qc_and_dublet_plot_location)
	for (j in 1:length(scRNA_set_names)){
		samples <- scRNA_sets[[j]]
		print(paste0(scRNA_data_path, samples[1], "/"))
		data_directory <- paste0(scRNA_data_path, samples[1], "/")
		sample1 <- make_seurat_object_and_doublet_removal(data_directory, samples[1])
		seu_list <- c()
		for (i in 2:length(samples)){
			data_directory <- paste0(scRNA_data_path, samples[i], "/")
			seu_list <- c(seu_list, make_seurat_object_and_doublet_removal(data_directory, samples[i]))
		}
		current_merge <- merge(sample1, y = seu_list, add.cell.ids = samples, project = scRNA_set_names[j])
		if (j==1){
			colon <- current_merge
		} else if (j>1){
			colon <- merge(colon, y = current_merge, project = "full_colon_project")
		}
	}
	setwd(analysis_parent_folder)
	colon[["percent.mt"]] <- PercentageFeatureSet(colon, pattern = "^MT-")
}
# 2) QC
if (2 %in% execute_steps){
	# create and set working directory to save qc plots
	if (!dir.exists(paste0(analysis_parent_folder, "all_samples_qc_plots"))){
		dir.create(paste0(analysis_parent_folder, "all_samples_qc_plots"))
	}
	setwd(paste0(analysis_parent_folder, "all_samples_qc_plots"))
	# make the standard seurat qc plots
	seurat_qc_plots(colon, sample_name)
	# Now subset the project (if not done already)
	colon <- subset(colon, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 5 & nCount_RNA < 10000)
	# leave the qc directory
	setwd(analysis_parent_folder)
}

metadata1 <- read.table(path_to_metadata, header = TRUE, sep = ",", stringsAsFactors=FALSE)
# 3) Add metadata
if (3 %in% execute_steps){
	metadata <- read.table(path_to_metadata, header = TRUE, sep = ",", stringsAsFactors=FALSE)
	# remove atac column
	metadata <- metadata[,colnames(metadata)[2:26]]
	colnames(metadata) <- c("Sample", colnames(metadata)[2:25])
	metadata <- metadata[metadata$Sample != "",]
	meta_data_types <- colnames(metadata)
	for (i in 2:length(meta_data_types)){
		identities <- colon[['orig.ident']]
		for (j in 1:length(metadata$Sample)){
			identities[identities==metadata$Sample[j]] <- metadata[j,meta_data_types[i]]
		}
		colon <- AddMetaData(colon, identities$orig.ident, col.name = meta_data_types[i])
	}
}
# 4) Normalize and scale data
if (4 %in% execute_steps){
	colon <- seurat_standard_normalize_and_scale(colon, TRUE, 1.0)
}
# 5) Plot UMAPs
if (5 %in% execute_steps){
	# Plot by clustering, sample, and disease state
	colon <- FindClusters(colon, resolution = 0.5)
	pdf(paste0("./UMAPclustering" , ".pdf"), onefile=F)
	print(DimPlot(colon, reduction = "umap", cols = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE)))
	dev.off()

	pdf(paste0("./UMAP_samples.pdf"), width = 12, onefile=F)
	print(DimPlot(colon, reduction = "umap", group.by = "orig.ident", cols = paletteDiscrete(values = unique(colon@meta.data$orig.ident), set = "stallion", reverse = FALSE)))
	dev.off()

	pdf(paste0("./UMAP_donor.pdf"), width = 12, onefile=F)
	print(DimPlot(colon, reduction = "umap", group.by = "Donor", cols = paletteDiscrete(values = unique(colon@meta.data$Donor), set = "stallion", reverse = FALSE)))
	dev.off()

	pdf(paste0("./UMAP_disease_state.pdf"), width = 6.5, onefile=F)
	print(DimPlot(colon, reduction = "umap", group.by = "DiseaseState",
		cols = c("#D51F26", "#D7CEC7", "#89288F", "#D7CEC7")) + theme_ArchR())
	dev.off()
	saveRDS(colon, "clustered_full_colon_proj_seurat.rds")
}

#==============================================================================================================================================
# annotation
dev.off()
rm(list=ls())
dir()
gc()
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ArchR)
library(viridis)
library(DoubletFinder)
library(tibble)
set.seed(1)
colon = readRDS("./clustered_full_colon_proj_seurat.rds")
DimPlot(colon, reduction = "umap",label = F,
        cols = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE))
DotPlot(colon, 
        features =c("EPCAM","CDH1","PTPRC","COL1A1","COL3A1","DES","MYL9","PECAM1","VWF"
                   ), 
        dot.scale = 8, 
        ) + RotatedAxis()
colon@meta.data$celltype2 = "NA"
Fibroblasts_clusters <- c(10,12,16,23,24)
Endo_clusters <- c(17)
immune_clusters <- c( 9,11,14,19,18,20,24)
epi_clusters <- c(0,1,2,3,4,5,6,7,8,13,15,21,22)
colon@meta.data[which(colon@meta.data$seurat_clusters %in% epi_clusters),'celltype2'] <- "Epithelial_cells"
colon@meta.data[which(colon@meta.data$seurat_clusters %in% Fibroblasts_clusters),'celltype2'] <- "Fibroblasts_cells"
colon@meta.data[which(colon@meta.data$seurat_clusters %in% Endo_clusters),'celltype2'] <- "Endothelial_cells"
colon@meta.data[which(colon@meta.data$seurat_clusters %in% immune_clusters),'celltype2'] <- "Immune_cells"
Fibroblasts <- DietSeurat(subset(colon, subset = seurat_clusters %in% Fibroblasts_clusters))
saveRDS(Fibroblasts, file = "colon_Fibroblasts_all_samples_initial.rds")
Endothelial <- DietSeurat(subset(colon, subset = seurat_clusters %in% Endo_clusters))
saveRDS(Endothelial, file = "colon_Endothelial_all_samples_initial.rds")
Idents(colon) = "celltype2"
DimPlot(colon, reduction = "umap",label = F,
        cols = paletteDiscrete(values = unique(colon@meta.data$celltype2), set = "stallion", reverse = FALSE))
cell.prop<-as.data.frame(prop.table(table(colon@meta.data$celltype2,colon$DiseaseState)))
colnames(cell.prop)<-c("cluster","sample","proportion")
p = ggplot(cell.prop,aes(sample,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  labs(x = '', y = 'Relative Abundance(%)') +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13)) +
  theme(legend.text = element_text(size = 11))
cell_type_cols <-paletteDiscrete(values = unique(colon@meta.data$celltype2), set = "stallion", reverse = FALSE)
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
FeaturePlot(colon,
            reduction = "umap",
            features = c("THBS2"),
            sort.cell = TRUE,
            pt.size = 1,
            label = F)
p3 = DotPlot(colon, 
        features =c("EPCAM","CDH1","COL1A1","COL3A1","DES","MYL9","PTPRC","PECAM1","EPAS1","VWF"
                   ), 
        dot.scale = 8, 
        ) + RotatedAxis()
p3
test1 <- subset(colon, DiseaseState == "Adenocarcinoma")
meta_data <- test1@meta.data %>% 
  rownames_to_column("Cell") %>% 
  select(Cell,celltype2) 
write.table(meta_data, 'meta_Adenocarcinoma.txt', sep='\t', quote=F, row.names=F)
write.table(as.matrix(test1@assays$RNA@data), 'count_Adenocarcinoma.txt', sep='\t', quote=F)
test2 <- subset(colon, DiseaseState == "Normal")
meta_data <- test2@meta.data %>% 
  rownames_to_column("Cell") %>% 
  select(Cell,celltype2) 
write.table(meta_data, 'meta_Normal.txt', sep='\t', quote=F, row.names=F)
write.table(as.matrix(test2@assays$RNA@data), 'count_Normal.txt', sep='\t', quote=F)
test3 <- subset(colon, DiseaseState == "Polyp")
meta_data <- test3@meta.data %>% 
  rownames_to_column("Cell") %>% 
  select(Cell,celltype2) 
write.table(meta_data, 'meta_Polyp.txt', sep='\t', quote=F, row.names=F)
write.table(as.matrix(test3@assays$RNA@data), 'count_Polyp.txt', sep='\t', quote=F)
test4 <- subset(colon, DiseaseState == "Unaffected")
meta_data <- test4@meta.data %>% 
  rownames_to_column("Cell") %>% 
  select(Cell,celltype2) 
write.table(meta_data, 'meta_Unaffected.txt', sep='\t', quote=F, row.names=F)
write.table(as.matrix(test4@assays$RNA@data), 'count_Unaffected.txt', sep='\t', quote=F)






