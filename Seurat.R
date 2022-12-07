rm(list = ls())
library(Seurat)
library(ggplot2)
library(psych)
library(magrittr)
library(dplyr)
library(patchwork)

data = read.csv('./merge_data.csv',header = T)
rownames(data) = data[,1]
data = data[,-1]
cui <- CreateSeuratObject(counts = data,project = "merge_data",min.cells = 3, min.features = 200)
table(cui@meta.data$orig.ident) 


cui <- subset(cui, subset = nFeature_RNA > 0)
dim(cui)
cui <- NormalizeData(cui, normalization.method = "LogNormalize", scale.factor = 10000)
cui <- FindVariableFeatures(cui, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(cui)
cui <- ScaleData(cui, features = all.genes)
cui <- RunPCA(cui,npcs = 50,features = VariableFeatures(object = cui),verbose = FALSE,approx=FALSE )

cui <- JackStraw(cui, num.replicate = 100)
cui <- ScoreJackStraw(cui, dims = 1:20)
ElbowPlot(cui)

cui <- FindNeighbors(cui, dims = 1:5)
cui <- FindClusters(cui, resolution = 0.1,random.seed = 256)
head(Idents(cui), 5)
table(cui@meta.data$seurat_clusters)

#UMAP
cui <- RunUMAP(cui, dims = 1:5)
head(cui@reductions$umap@cell.embeddings) 
pdf("./celltype_umap.pdf")
DimPlot(cui, reduction = "umap",label = T)
dev.off()


#T-SNE
cui <- RunTSNE(cui,dims = 1:5,check_duplicates = FALSE)
head(cui@reductions$tsne@cell.embeddings)
pdf("./celltype_tsne.pdf")
DimPlot(cui, reduction = "tsne",label = T)
dev.off()


cui.markers <- FindAllMarkers(cui, only.pos = FALSE,test.use = "wilcox")
top3 <- cui.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
pdf("./cellmarker_vlnplot.pdf")
VlnPlot(cui, features = top3$gene,pt.size=0)
dev.off()
pdf("./cellmarker_tsne.pdf")
FeaturePlot(cui,reduction = "tsne",features = top3$gene)
dev.off()

library(SingleR)
head(cui@meta.data)
str(cui)
cui <- RunTSNE(cui, dims = 1:5,check_duplicates = FALSE)
cui <- RunUMAP(cui, dims = 1:5)
load(file="./hpca.se.RData")
head(colnames(hpca.se))
head(rownames(hpca.se))
testdata <- GetAssayData(cui, slot="data")
dim(testdata)
testdata[1:30,1:4]
clusters <- cui@meta.data$seurat_clusters
table(clusters)
cellpred <- SingleR(test = testdata,  
                    ref = hpca.se, 
                    labels = hpca.se$label.main,
                    method = "cluster", 
                    clusters = clusters)

celltype = data.frame(ClusterID = rownames(cellpred),
                      celltype = cellpred$labels,
                      stringsAsFactors = F)


new.cluster.ids <- cellpred$pruned.labels
names(new.cluster.ids) <- levels(cui)
cui <- RenameIdents(cui,new.cluster.ids)
cui[['cell_type']]<-cui@active.ident
pdf("./celltype_anno_tsne.pdf")
DimPlot(cui, reduction = "tsne",label = T)
dev.off()
pdf("./celltype_anno_umap.pdf")
DimPlot(cui, reduction = "umap",label = T)
dev.off()


save(cui,file = './seurat_singleR_mergedata.Rdata')


