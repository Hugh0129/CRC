# process data
import os, sys 
os.getcwd()
os.listdir(os.getcwd())
import loompy as lp
import numpy as np
import scanpy as sc
x=sc.read_csv("scRNA_epi.csv")
row_attrs = {"Gene": np.array(x.var_names),}
col_attrs = {"CellID": np.array(x.obs_names)}
lp.create("sample.loom",x.X.transpose(),row_attrs,col_attrs)
exit   

#=========================================================================================
# pyscenic
pyscenic grn \
--num_workers 20 \
--output adj.sample.tsv \
--method grnboost2 \
sample.loom \
hs_hgnc_tfs.txt 


pyscenic ctx \
adj.sample.tsv \
hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
--annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname sample.loom \
--mode "dask_multiprocessing" \
--output reg.csv \
--num_workers 6 \
--mask_dropouts


pyscenic aucell \
sample.loom \
reg.csv \
--output sample_SCENIC.loom \
--num_workers 6

#=====================================================================================================
dev.off()
rm(list=ls())
dir()
gc()
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(ggplot2)
library(pheatmap)
library(Seurat)
library(ggheatmap)
library(reshape2)
sce_SCENIC <- open_loom('../data/sample_SCENIC.loom')         
regulons_incidMat <- get_regulons(sce_SCENIC, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
class(regulons)
regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(sce_SCENIC)
human_data <- readRDS("CAF.RDS")
cellinfo <- human_data@meta.data[,c('spefic_celltype','DiseaseState',"nFeature_RNA","nCount_RNA")]
colnames(cellinfo)=c('celltype', 'group','nGene' ,'nUMI')
cellTypes <-  as.data.frame(subset(cellinfo,select = 'celltype'))
selectedResolution <- "celltype"
sub_regulonAUC <- regulonAUC
rss <- calcRSS(AUC=getAUC(sub_regulonAUC),
               cellAnnotation=cellTypes[colnames(sub_regulonAUC),
                                        selectedResolution])
rss=na.omit(rss)
ex = plotRSS(
  rss,
  labelsToDiscard = NULL,
  zThreshold = 2,
  cluster_columns = FALSE,
  order_rows = TRUE,
  thr = 0.1,
  varName = "cellType",
  col.low = "grey90",
  col.mid = "darkolivegreen3",
  col.high = "darkgreen",
  revCol = FALSE,
  verbose = TRUE
)
ex
rss_data <- rssPlot$plot$data
rss_data<-dcast(rss_data, 
                Topic~rss_data$celltype,
                value.var = 'Z')
rownames(rss_data) <- rss_data[,1]
rss_data <- rss_data[,-1]
colnames(rss_data)
col_ann <- data.frame(group= c(rep("Unaffected_Polyp_CAF_infla_state1",1),
                               rep("Pre_CAF_infla_state1",1),
                               rep("NF_state1",1),
                               rep("THBS2_CAF_infla_state1",1)))
rownames(col_ann) <- colnames(rss_data)
groupcol <- c("#D9534F", "#96CEB4", "#CBE86B", "#EDE574")
names(groupcol) <- c("Unaffected_Polyp_CAF_infla_state1",
                     "Pre_CAF_infla_state1",
                     "NF_state1", 
                     "THBS2_CAF_infla_state1")
col <- list(group=groupcol)

text_columns <- sample(colnames(rss_data),0)
rownames(rss_data)
p <- ggheatmap(rss_data,color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
               cluster_rows = T,cluster_cols = F,scale = "row",
               annotation_cols = col_ann,
               annotation_color = col,
               legendName="Relative value",
               text_show_cols = text_columns)
p
























































































