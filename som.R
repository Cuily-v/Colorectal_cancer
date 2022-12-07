rm(list = ls())
library(kohonen)
library(Seurat)

load(file='./seurat_singleR_mergedata.Rdata')
counts_matrix = as.matrix(GetAssayData(cui, slot="counts"))
CNV_sample = as.matrix(read.csv(file = './CNV_sample.csv'))

counts_matrix = t(counts_matrix )
counts_matrix = counts_matrix[CNV_sample,]
dim(counts_matrix)
data = t(counts_matrix)

data <- as.matrix(t(scale(t(data)))) 
head(data) 
N = dim(data)[2]
size = round(sqrt(5*sqrt(N)))


som_grid <- somgrid(xdim = size, 
                    ydim = size, 
                    topo="rectangular",
                    neighbourhood.fct ="gaussian") 

weights <- matrix(rnorm((dim(data)[1])*(dim(data)[2])),nrow = dim(data)[1],ncol = dim(data)[2])
dim(weights)

som_model <- som(data,                 
                 grid=som_grid,                      
                 rlen = 100,                     
                 #alpha = 0.5,                      
                 #radius =  3,                     
                 keep.data = TRUE,                      
                 user.weights = weights )
plot(som_model, type = "changes")



coolBlueHotRed <- function(n, alpha = 0.5) {
  rainbow(n, end = 4/6, alpha=alpha)[n:1]
}
plot(som_model, type = "counts", main="Node Counts", palette.name=coolBlueHotRed)
plot(som_model, type = "quality", 
     main="Node Quality/Distance", 
     palette.name=coolBlueHotRed)
plot(som_model, type = "mapping",
     border = "gray",
     shape = "straight",
     palette.name=coolBlueHotRed)
plot(som_model, type="dist.neighbours",
     main = "SOM neighbour distances", 
     palette.name=grey.colors)


table(som_model$unit.classif)
som_model_code_class = data.frame(name=rownames(data), code_class=som_model$unit.classif)
head(som_model_code_class)

mydata <- as.matrix(as.data.frame(som_model$codes)) 
dim(mydata)
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata, centers=i)$withinss)
par(mar=c(5.1,4.1,4.1,2.1))
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares", main="Within cluster sum of squares (WCSS)")


# Form clusters on grid
# use hierarchical clustering to cluster the codebook vectors
som_cluster <- cutree(hclust(dist(mydata)), 5)
# Colour palette definition
cluster_palette <- function(x, alpha = 0.5) {
  n = length(unique(x)) * 2
  rainbow(n, start=2/6, end=6/6, alpha=alpha)[seq(n,0,-2)]
}

cluster_palette_init = cluster_palette(som_cluster)
bgcol = cluster_palette_init[som_cluster]
#show the same plot with the codes instead of just colours
plot(som_model, type="codes", bgcol = bgcol, main = "Clusters", codeRendering="lines")
add.cluster.boundaries(som_model, som_cluster)
som_model_code_class_cluster = som_model_code_class
som_model_code_class_cluster$cluster = som_cluster[som_model_code_class$code_class]
head(som_model_code_class_cluster)

