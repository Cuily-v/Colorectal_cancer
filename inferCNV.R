# cuilingyu
# 2022-7-28


rm(list = ls())
library(rjags)
library(infercnv)
library(dplyr)
library(Seurat)
library(SingleR)
library(magrittr)

load(file='./seurat_singleR_mergedata.Rdata')
counts_matrix = as.matrix(GetAssayData(cui, slot="counts"))

ann=as.matrix(cui@active.ident)
ann=cbind(rownames(ann),ann)
write.table(ann, './ann.txt',  quote = FALSE,row.names = FALSE,col.names = FALSE,sep='\t')

refgene <- read.table("./refGene.txt", sep = "\t")
refgene<-refgene[,c(3,5,6,13)]
refgene <- refgene%>% distinct(V13, .keep_all = T)
rownames(refgene)<-refgene[,4]
names=as.matrix(intersect(row.names(counts_matrix),rownames(refgene)))
gene_infor=na.omit(refgene[names,])
gene_infor <- gene_infor[,c(4,1,2,3)]
write.table(gene_infor, './gene_infor.txt',  quote = FALSE,row.names = FALSE,col.names = FALSE,sep='\t')

counts_matrix=na.omit(counts_matrix[names,])
dim(counts_matrix)

obj = CreateInfercnvObject(
  raw_counts_matrix=counts_matrix,
  annotations_file='./ann.txt',
  ref_group_names=NULL,
  gene_order_file='./gene_infor.txt',
  delim="\t")

infercnv_obj = infercnv::run(obj,
                             cutoff = 0.1, 
                             out_dir = './NULL_ref',
                             cluster_by_groups = FALSE, 
                             denoise = TRUE,
                             output_format = "pdf",
                             HMM = FALSE)


