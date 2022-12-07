rm(list=ls())
library(openxlsx)

scGene_exp <- read.csv(file = './gene_exp.csv',header = T)[,-1]
scATAC_seq <- read.csv(file = './scATAC-seq.csv',header = T)
same_gene_names  =  as.data.frame(intersect(scGene_exp$rowname,scATAC_seq$X))
colnames(same_gene_names) = c("X") 


scATAC_seq = merge(scATAC_seq, same_gene_names,by = "X")

library(dplyr) 
scGene_exp = scGene_exp %>% group_by(rowname) %>% summarise_each(funs(mean))
colnames(same_gene_names) = c("rowname") 
scGene_exp = merge(scGene_exp, same_gene_names,by = "rowname")

rownames(scATAC_seq)=scATAC_seq [,1]
scATAC_seq = scATAC_seq [,-1]
rownames(scGene_exp)=scGene_exp[,1]
scGene_exp= scGene_exp[,-1]
scGene_exp = round(scGene_exp,2)
scATAC_seq  = round(scATAC_seq ,2)

data = cbind(scATAC_seq,scGene_exp)
write.csv(data, file = './merge_data.csv', row.names = T)

