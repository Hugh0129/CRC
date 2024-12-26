dev.off()
rm(list=ls())
dir()
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE) 
library(Mfuzz)
library(ggplot2)
library(plyr)
library(ggord)
library(yyplot)
library(limma)
library(stats)
library(ggsci)
library(ggThemeAssist)
library(RColorBrewer)
expr_df = read.csv('Expr_tpm_ok.csv', header = T,row.names=1, sep = ',', quote = '', stringsAsFactors = FALSE, check.names = FALSE)
meta_df = read.csv('meta_data.csv',row.names=1,header = T,sep = ',',quote = '', check.names = FALSE)
expr_df[1:3,1:4]
Expr = as.data.frame(expr_df)
l = read.table("gene_list.txt",header = F)
Expr = Expr[l$V1,]
expr_df = as.data.frame(t(Expr))
expr_df$label = substr(rownames(expr_df),0,1)
expr_df <- expr_df %>% dplyr::select('label', everything())
expr_df[1:3,1:4]
df2<-aggregate(expr_df[,colnames(expr_df)[2:ncol(expr_df)]],by=list(expr_df$label),mean,na.rm= TRUE)
row.names(df2)<-df2[,1]
df3<-data.frame(t(df2[,-1]))
df3 = df3 %>% dplyr::select('N','P','C')
## 2.1 Filtering----
df3a<-as.matrix(df3)
df3Ex<- ExpressionSet(assayData = df3a)
df3F <- filter.NA(df3Ex,thres = 0.25)
df3F <- fill.NA(df3F,mode = 'mean')
df3F <- filter.std(df3F, min.std = 0)
## 2.2 Standardisation----
df3F <- standardise(df3F)
## 2.3 Setting of parameters for FCM clustering----
set.seed(1)
cluster_num = 6
mestimate(df3F)
cl <- mfuzz(df3F,c=cluster_num,m=mestimate(df3F))
mfuzz.plot2(df3F, cl=cl,mfrow=c(3,2),centre=TRUE,x11=F,centre.lwd=0.2, 
            time.labels = colnames(df3))
dir.create(path="mfuzz",recursive = TRUE)
for(i in 1:cluster_num){
  potname<-names(cl$cluster[unname(cl$cluster)==i])
  write.csv(cl[[4]][potname,i],paste0("mfuzz","/mfuzz_",i,".csv"))
}

cluster_size <- cl$size
names(cluster_size) <- 1:cluster_num
cluster_size
protein_cluster <- cl$cluster
protein_cluster <- cbind(df3[names(protein_cluster), ], protein_cluster)
head(protein_cluster)
write.table(protein_cluster, 'protein_cluster.txt', sep = '\t', col.names = NA, quote = FALSE)



