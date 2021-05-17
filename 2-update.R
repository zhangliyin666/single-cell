#1.Preparation for data
library(dplyr)
library(Seurat)
library(cowplot)
library(patchwork)
library(ggplot2)
#put all data into one dictionary
samples=list.files("/home/data/vip9t07/lungfibro/adult/data/")
samples
#2.Create Seurat object
#loop for reading data and create seurat object
sceList = lapply(samples,function(pro){
  folder=file.path("/home/data/vip9t07/lungfibro/adult/data/",pro)
  CreateSeuratObject(counts = Read10X(folder),
                     project = pro )
})
#manually add groups
sceList[[1]]@meta.data$stim = "fibrosis"
sceList[[2]]@meta.data$stim = "fibrosis"
sceList[[3]]@meta.data$stim = "control"
sceList[[4]]@meta.data$stim = "control"

#merge the objects
sce.big <- merge(sceList[[1]],
                 y = c(sceList[[2]],sceList[[3]],sceList[[4]]),
                 add.cell.ids = samples,
                 project = "ls_12")
sce.big

table(sce.big$orig.ident)

#save(sce.big,file = 'sce.big.merge.ls_12.Rdata')
#3.Filteration and quantification control
raw_sce=sce.big
#mito gene
rownames(raw_sce)[grepl('^MT-',rownames(raw_sce),ignore.case = T)]
#ribo gene
rownames(raw_sce)[grepl('^RP[SL]',rownames(raw_sce),ignore.case = T)]
#calculate the proportion of mito and ribo gene
raw_sce[["percent.mt"]] <- PercentageFeatureSet(raw_sce, pattern = "^MT-")

rb.genes <- rownames(raw_sce)[grep("^RP[SL]",rownames(raw_sce))]
C<-GetAssayData(object = raw_sce, slot = "counts")
percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
raw_sce <- AddMetaData(raw_sce, percent.ribo, col.name = "percent.ribo")

plot1 <- FeatureScatter(raw_sce, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "stim")
plot2 <- FeatureScatter(raw_sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

VlnPlot(raw_sce, features = c("percent.ribo", "percent.mt"), ncol = 2,pt.size = 0)
VlnPlot(raw_sce, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2,pt.size = 0)
VlnPlot(raw_sce, features = c("percent.ribo", "nCount_RNA"), ncol = 2,pt.size = 0)

raw_sce 
raw_sce1 <- subset(raw_sce, subset = nFeature_RNA > 200 & nCount_RNA > 1000 & percent.mt < 20)
raw_sce1

sce=raw_sce1
#perform standard preprocessing (log-normalization), and identify variable features individually for each
pancreas.list <- SplitObject(sce, split.by = "stim")
for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000, 
                                             verbose = FALSE)
}
#4.Dimension reduction and clustering
#identify anchors using the function, which takes a list of Seurat objects as input
pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, dims = 1:30)#10~50 are OK
#save(pancreas.anchors,file = 'pancreas.anchors.Rdata')
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
# switch to integrated assay. The variable features of this assay are automatically set during
# IntegrateData
DefaultAssay(pancreas.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
pancreas.integrated <- FindNeighbors(pancreas.integrated, reduction = "pca", dims = 1:20)
pancreas.integrated <- FindClusters(pancreas.integrated, resolution = 0.5)
DimPlot(pancreas.integrated,label = 1)
DimPlot(pancreas.integrated, reduction = "umap", split.by = "stim")
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "stim")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "orig.ident",label = TRUE, repel = TRUE) 
p1 + p2

#save(pancreas.integrated,file = 'pancreas.integrated.Rdata')
#5.Cluster classification
DefaultAssay(pancreas.integrated) <- "RNA"


#6.Subset interested clusters
Idents(pancreas.integrated)
levels(pancreas.integrated)
head(pancreas.integrated@meta.data)

genes_to_check = c('PECAM1','VWF','CDH5')

DotPlot(pancreas.integrated, # group.by = 'seurat_clusters',
        features = unique(genes_to_check)) + RotatedAxis()

p1=DimPlot(pancreas.integrated, reduction = 'umap', group.by = 'seurat_clusters',
           label = TRUE, pt.size = 0.5) + NoLegend()
p2=DotPlot(pancreas.integrated, group.by = 'seurat_clusters',
           features = unique(genes_to_check)) + RotatedAxis()
library(patchwork)
p1+p2

sce = pancreas.integrated
endo1 = sce[,sce@meta.data$seurat_clusters %in% c(2,9,16)]

endo1 <- NormalizeData(endo1, normalization.method = "LogNormalize", scale.factor = 1e4) 
endo1 <- FindVariableFeatures(endo1, selection.method = 'vst', nfeatures = 2000)
endo1 <- ScaleData(endo1, vars.to.regress = "percent.mt")
endo1 <- RunPCA(endo1, features = VariableFeatures(object = endo1)) 

endo1 <- FindNeighbors(endo1, dims = 1:10)
endo1 <- FindClusters(endo1, resolution = 1 )
# Look at cluster IDs of the first 5 cells
head(Idents(endo1), 5)
table(endo1$seurat_clusters) 
endo1 <- RunUMAP(endo1, dims = 1:10)
DimPlot(endo1, reduction = 'umap')

genes_to_check = c('PECAM1','VWF','CDH5')
DotPlot(endo1, group.by = 'seurat_clusters',
        features = unique(genes_to_check)) + RotatedAxis()

p1=DimPlot(endo1, reduction = 'umap', group.by = 'seurat_clusters',
           label = TRUE, pt.size = 0.5) + NoLegend()
p2=DotPlot(endo1, group.by = 'seurat_clusters',
           features = unique(genes_to_check)) + RotatedAxis()
p1+p2

DimPlot(endo1,label = 1)
DimPlot(endo1, reduction = "umap", split.by = "stim",label = 1)
p1 <- DimPlot(endo1, reduction = "umap", group.by = "stim")
p2 <- DimPlot(endo1, reduction = "umap", group.by = "orig.ident",label = TRUE, repel = TRUE) 
p1 + p2

plots <- VlnPlot(endo1, features = c('RHOJ'), split.by = "stim", group.by = "seurat_clusters", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

#save(endo1,file = 'endo1.Rdata')
