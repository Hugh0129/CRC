dev.off()
rm(list=ls())
dir()
gc()
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE) 
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
library(Ktplots)
CAF <- readRDS("./CAF_need.RDS")
#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(CAF@assays$RNA@counts), 'sparseMatrix')
pd <- CAF@meta.data
pd$cell.type<-Idents(CAF)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
cds <- new_cell_data_set(expression_data = data,
                         cell_metadata = pd,
                         gene_metadata = fData)
#NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 100)  
plot_pc_variance_explained(cds)   
cds <- align_cds(cds, alignment_group = "SampleNameOnly")
cds <- reduce_dimension(cds,preprocess_method = "PCA") 
cds <- reduce_dimension(cds, reduction_method="tSNE")
cds <- cluster_cells(cds)   
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(CAF, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed 
plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters") 
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "spefic_celltype",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           group_label_size=4,
           cell_size=1.5)          
cds = order_cells(cds)  
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1.5)


