rm(list=ls())
library(infercnv)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")

infercnv_obj = readRDS("./NULL_ref/run.final.infercnv_obj")
expr <- infercnv_obj@expr.data
test_loc <- infercnv_obj@observation_grouped_cell_indices
test_loc <- test_loc$Epithelial_cells
anno.df=data.frame(
  CB=c(colnames(expr)[test_loc]),
  class=c(rep("test",length(test_loc)))
)
head(anno.df)


set.seed(256)
kmeans.result <- kmeans(t(expr), 5)
kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
kmeans_df$CB=rownames(kmeans_df)
kmeans_df=kmeans_df%>%inner_join(anno.df,by="CB") 
kmeans_df_s=arrange(kmeans_df,kmeans_class)
rownames(kmeans_df_s)=kmeans_df_s$CB
kmeans_df_s$CB=NULL
kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class)
head(kmeans_df_s)
color_v=RColorBrewer::brewer.pal(8, "Dark2")[1:5]
names(color_v)=as.character(1:5)


expr.scale <- scale(t(expr))
tmp1 <- sweep(expr.scale, 2, apply(expr.scale, 2, min),'-') 
tmp2 <- apply(expr.scale, 2, max) - apply(expr.scale,2,min) 

expr_1 <- t(2*sweep(tmp1, 2, tmp2, "/")-1)  
CNV_score <- as.data.frame(colSums(expr_1 * expr_1)) 
colnames(CNV_score)="CNV_score"
CNV_score <- rownames_to_column(CNV_score, var='cell')


CNV_score$CB=CNV_score$cell
kmeans_df_s$CB=rownames(kmeans_df_s)
CNV_score=CNV_score%>%inner_join(kmeans_df_s,by="CB")
CNV_score%>%ggplot(aes(kmeans_class,CNV_score))+geom_violin(aes(fill=kmeans_class),color="NA")+
  scale_fill_manual(values = color_v)+
  geom_boxplot(width = 0.1, outlier.size = 0.5 ,notch = T,show.legend=T)+
  geom_hline(aes(yintercept=500),linetype=500,col="red")+
  theme_bw()
ggsave("./CNV_5_class.pdf",width = 10,height = 6,units = "cm")


rownames(CNV_score) = CNV_score[,1]
CNV_sample = CNV_score[CNV_score$CNV_score > 500,][,1]
write.csv(CNV_sample,file = './CNV_sample.csv',row.names = F)

