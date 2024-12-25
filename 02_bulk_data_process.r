dev.off()
rm(list=ls())
dir()
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
library(ggplot2)
library(plyr)
library(ggord)
library(yyplot)
library(org.Hs.eg.db)
library(clusterProfiler)
library(limma)
expr_df_tpm = read.csv('ensembl_matrix_TPM.csv', header = T,row.names=1, sep = ',', quote = '', stringsAsFactors = FALSE, check.names = FALSE) 
expr_df_count = read.csv('ensembl_matrix_TPM_count.csv', header = T,row.names=1, sep = ',', quote = '', stringsAsFactors = FALSE, check.names = FALSE) 
symbol_tpm = as.character(substr(rownames(expr_df_tpm),1,15))
eg1_tpm = bitr(symbol_tpm, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
expr_df_tpm$ENSEMBL =  substr(rownames(symbol_tpm),1,15)
Expr_tpm = merge(eg1_tpm, expr_df_tpm, by="ENSEMBL")
Expr_tpm$ENSEMBL = NULL
Expr_tpm = avereps(Expr_tpm[,-1],ID = Expr_tpm$SYMBOL) 
Expr_tpm = as.data.frame(Expr_tpm)
symbol_count = as.character(substr(rownames(expr_df_count),1,15))
eg1_count = bitr(symbol_count, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
expr_df_count$ENSEMBL = substr(rownames(expr_df_count),1,15)
Expr_count = merge(eg1_count, expr_df_count, by="ENSEMBL")
Expr_count$ENSEMBL = NULL
Expr_count = avereps(Expr_count[,-1],ID = Expr_count$SYMBOL)
Expr_count = as.data.frame(Expr_count)
Expr_count[1:3,1:4]
Expr_count = Expr_count[rowMeans(Expr_count)>10,] 
table(rownames(Expr_tpm) %in% rownames(Expr_count))
Expr_tpm_ok = Expr_tpm[(rownames(Expr_tpm) %in% rownames(Expr_count)),]
Expr_count_ok = Expr_count[rownames(Expr_count)%in%rownames(Expr_tpm_ok),]
table(rownames(Expr_count_ok) == rownames(Expr_tpm_ok))
write.csv(Expr_tpm_ok, file = "Expr_tpm_ok.csv", quote=F, row.names=T)
write.csv(Expr_count_ok, file = "Expr_count_ok.csv", quote=F, row.names=T)











