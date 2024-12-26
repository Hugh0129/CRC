dev.off()
rm(list=ls())
dir()
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE) 
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
data = read.csv('Expr_tpm_ok.csv', header = T,row.names=1, sep = ',', quote = '', stringsAsFactors = FALSE, check.names = FALSE) 
data = as.data.frame(t(data))
pdata = read.csv('meta_data.csv', header = T,row.names=1, sep = ',', quote = '', stringsAsFactors = FALSE, check.names = FALSE)  # quote去除引号
data = data[rownames(data)%in%rownames(pdata),]
pdata = pdata[rownames(pdata)%in%rownames(data),]
data = data[rownames(pdata),]
all(rownames(data)==rownames(pdata))
data_need = data[,c("THBS2","CD36")]
table(plot_data$type)
plot_data = cbind(pdata,data_need)
plot_data = plot_data[,c(5,7)]
colnames(plot_data) = c("group","Retive_Abundance")
levels(factor(plot_data$group))
plot_data$group <- factor(plot_data$group,levels = c("normal","adenoma","adenocarcinoma"))
p<- ggplot(data=plot_data)+ 
  geom_boxplot(mapping=aes(x=group,y=Retive_Abundance,colour = group ),
               alpha = 0.5,
               size=0.75,
               width = 0.6)+ 
  geom_jitter(mapping=aes(x=group,y=Retive_Abundance,colour = group), 
              alpha = 1,size=1)+
  scale_color_manual(limits=c("normal","adenoma","adenocarcinoma"), 
                     values=c("#85B22E","#5F80B4","#922927"))+ 
  geom_signif(mapping=aes(x=group,y=Retive_Abundance),
              comparisons = list(c("normal", "adenoma"), 
                                 c("adenoma", "adenocarcinoma"),
                                 c("normal", "adenocarcinoma")),
              map_signif_level=T, 
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), 
              y_position = c(21,22,23), 
              size=0.75, 
              textsize = 4, 
              test = wilcox.test)+ 
  theme_classic(  
    base_line_size = 0.75
  )+
  labs(title="CD36",x="",y="Gene expression(TPM)")+ 
  theme(plot.title = element_text(size = 15,
                                  colour = "black",
                                  hjust = 0.5),
        axis.title.y = element_text(size = 15, 
                                    color = "black",
                                    face = "plain", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_text(color="black", 
                                    size=15, 
                                    face="plain"),
        legend.text = element_text(color="black", 
                                   size = 10, 
                                   face = "plain"),
        axis.text.x = element_text(size = 13, 
                                   color = "black",
                                   face = "plain",
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0), 
        axis.text.y = element_text(size = 13,  
                                   color = "black",
                                   face = "plain", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  )
p= p +  ylim(0,20)#+ labs(title="Volcanoplot",x=expression(logFC),y=expression(-log10(P.Value)))
p









