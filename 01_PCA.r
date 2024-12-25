

dev.off()
rm(list=ls())
dir()

getwd()
setwd("G:/work/who/me/zheyiyangben/00_pca/results/")
getwd()
dir()

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor-0

library(ggplot2)
library(plyr)
# devtools::install_github('fawda123/ggord')
library(ggord)
# devtools::install_github("GuangchuangYu/yyplot")
library(yyplot)


expr_df = read.csv('G:/work/who/me/zheyiyangben/00_count2tpm/results/ensembl_matrix_TPM.csv', header = T,row.names=1, sep = ',', quote = '', 
             stringsAsFactors = FALSE, check.names = FALSE)  # quote去除引号
meta_df = read.csv('G:/work/who/me/zheyiyangben/raw_data/meta_data.csv',row.names=1,
                 header = T,sep = ',',quote = '', check.names = FALSE)


#查看前3个基因在前4个sample中的表达矩阵
expr_df[1:3,1:4]
expr_df = expr_df[rowMeans(expr_df)>10,] # 自定义



symbol = as.character(rownames(expr_df))
eg1 = bitr(symbol, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")

expr_df$ENSEMBL = rownames(expr_df)
Expr = merge(eg1, expr_df, by="ENSEMBL")

Expr$ENSEMBL = NULL
#对重复基因名取平均表达量，然后将基因名作为行名
Expr = avereps(Expr[,-1],ID = Expr$SYMBOL) # 自定义

Expr = as.data.frame(Expr)



expr_df = as.data.frame(t(Expr))

#查看样本信息前3行
head(meta_df, n=3)


# 做PCA
#用`prcomp`进行PCA分析
pca.results <- prcomp(expr_df, center = TRUE, scale. = TRUE)

#定义足够多的颜色，用于展示分组
mycol <- c("#223D6C","#D20A13","#088247","#FFD121","#11AA4D","#58CDD9","#7A142C",
           "#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")

#用ggord画基本PCA图
ggord(pca.results, grp_in = meta_df$group, repel=TRUE,
      xlims = c(-110,110),
      ylims = c(-60,110),
      ellipse = FALSE, #不显示置信区间背景色
      size = 2, #样本的点大小
      alpha=0.5, #设置点为半透明，出现叠加的效果
      #如果用自定义的颜色，就运行下面这行
      cols = mycol[1:length(unique(meta_df$group))],
      arrow = NULL,txt = NULL) + #不画箭头和箭头上的文字
  theme(panel.grid =element_blank()) # + #去除网格线
  
  #用yyplot添加置信区间圆圈
  # geom_ord_ellipse(ellipse_pro = .95, #设置置信区间
  #                  size=1.5, #线的粗细
  #                  lty=1 ) #实线




