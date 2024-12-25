dev.off()
rm(list=ls())
dir()
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE) 
library(limma)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(data.table)
library(clusterProfiler)
library(DOSE)
library(GO.db)
library(org.Hs.eg.db)
library(topGO)
library(GSEABase)
library(stringr)
a = read.csv('count_ok.csv', header = T,row.names=1, sep = ',', quote = '', stringsAsFactors = FALSE, check.names = FALSE) 
a=as.data.frame(a)
a[1:4,1:4]
Expr = a[,c(1:10,21:30)]
head(Expr)
pdata = read.csv('meta_data.csv',row.names=1,header = T,sep = ',',quote = '',check.names = FALSE)
pdata <- as.data.frame(pdata)
pdata[1:3,]
pdata = pdata[colnames(Expr),]
aa = 0
ak <- function(x) {
  if (as.numeric(x) > aa) {
    x = as.factor("high")
    return(x)
  } else {
    x = as.factor("low")
    return(x)
  }
}
pp <- data.frame(sapply(pdata$type, function(x) ak(x)))
rownames(pp) = rownames(pdata)
colnames(pp) = "type" 
#Step1
#limma
group_list = as.character(pp[, 1])
table(group_list)
# Multidimensional scaling (MDS) plot
plotMDS(Expr, col = as.numeric(group_list))
# Step2
design <- model.matrix(~0+factor(group_list))
levels(factor(group_list))
colnames(design) = levels(factor(group_list))
rownames(design) = rownames(pdata)
# Step3
DGElist = DGEList(counts = Expr, group = group_list)
keep_gene = rowSums( cpm(DGElist) > 1 ) >= 10
table(keep_gene)
DGElist = DGElist[keep_gene, , keep.lib.sizes = FALSE]
# Step4
DGElist = calcNormFactors(DGElist, method = 'TMM') ##TMM
v = voom(DGElist, design, plot = TRUE, normalize = "quantile") 
fit = lmFit(v, design)
cont.matrix = makeContrasts(contrasts = c('high-low'), levels = design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)
plotSA(fit2, main="Final model: Mean-variance trend")
# Step5
nrDEG_limma_voom = topTable(fit2, coef = 'high-low', n = Inf, adjust="BH")
nrDEG_limma_voom = na.omit(nrDEG_limma_voom)
head(nrDEG_limma_voom)
write.csv(nrDEG_limma_voom, file = "mRNA_nrDEG_matrix_CvsN.csv", quote=F, row.names=T)
q = 0.05 
LogFoldChange= 1 
nrDEG_limma_voom_signif = nrDEG_limma_voom[(nrDEG_limma_voom$adj.P.Val < q & 
                                              (nrDEG_limma_voom$logFC>LogFoldChange | 
                                                 nrDEG_limma_voom$logFC<(-LogFoldChange))),]
dim(nrDEG_limma_voom_signif)
nrDEG_limma_voom_signif = nrDEG_limma_voom_signif[order(nrDEG_limma_voom_signif$logFC),]
write.csv(nrDEG_limma_voom_signif, file="mRNA_nrDEG_signif.csv", quote=F, row.names=T)
id = rownames(nrDEG_limma_voom_signif)
matrix_diff =  Expr[id,]
write.csv(matrix_diff, file="matrix_mRNA_diff.csv", quote=F, row.names=T)  


















