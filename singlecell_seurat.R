library(dplyr)
library(Seurat)
library(cowplot)

mj1 <- Read10X(data.dir = '1-MJ-PH4D-WT/')
mj2 <- Read10X(data.dir = '2-MJ-PH4D-KO/')
cyt1 <- Read10X(data.dir = '3-CYT -PH0D-WT/')
#cyt2 <- Read10X(data.dir = '4-CYT- PH4D-WT/')

object1 <- CreateSeuratObject(counts = mj1, project = "pbmc3k", min.cells = 3, min.features = 200)
object2 <- CreateSeuratObject(counts = mj2, project = "pbmc3k", min.cells = 3, min.features = 200)
object3 <- CreateSeuratObject(counts = cyt1, project = "pbmc3k", min.cells = 3, min.features = 200)
#object2 <- CreateSeuratObject(counts = cyt2, project = "pbmc3k", min.cells = 3, min.features = 200)

object1@meta.data$stim='ph4d-wt'  #MJ-PH4D-WT
object2@meta.data$stim='ph4d-ko'  #MJ-PH4D-KO
object3@meta.data$stim='ph0d-wt'  #CYT -PH0D-WT
#object2@meta.data$stim='ph4d-wt'  #CYT -PH4D-WT

data = merge(object1,c(object2,object3))
data[['percent.mt']] <- PercentageFeatureSet(data,pattern='^mt-')
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
table(data$stim)
ifnb = data
ifnb.list <- SplitObject(ifnb, split.by = "stim")
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
#FindAnchors、降维
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = 1)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = 1)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20,min.dist = 0.4)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.8)
DimPlot(immune.combined,label = 1)

immune.combined <- subset(immune.combined,idents = c(0:17))
uncertain <- subset(immune.combined,idents = c(7,9))
DimPlot(uncertain,label = 1)

p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
plot_grid(p1)
DimPlot(immune.combined, reduction = "umap", split.by = "stim")


saveRDS(immune.combined, file = "immunecombined.rds")

DefaultAssay(immune.combined) <- "RNA"
VlnPlot(immune.combined,features = c('Adm','Top2a','Gja4','Cd36','Ednrb','Edn1','Vwf','Pdpn'),pt.size = 0)#
VlnPlot(immune.combined,features = c('Adm','Pgf','Cxcr4','Esm1','Apln','Angpt2'),pt.size = 0) #tip cell
VlnPlot(immune.combined,features = c('Pxdn','Plxnd1','Lamb1','Lama4','Nid2','Hspg2','Igfbp3'),pt.size = 0) #breach cell
VlnPlot(immune.combined,features = c('Ptprc','Pecam1','Mki67','Cd200'),pt.size = 0)
VlnPlot(immune.combined,features = c('Vwf','Vcam1','Nr2f2','Nrp2','Ackr3'),pt.size = 0) #vein
VlnPlot(immune.combined,features = c('Gja4','Sox17','Cxcl12','Efnb2','Sema3c','Gja5'),pt.size = 0) #artery
VlnPlot(immune.combined,features = c('Alb','Afp','Albumin','Acly','Apoa','Asl','Ass'),pt.size = 0) #hepatocyte
VlnPlot(immune.combined, features = 'Adm',pt.size = 0)#tip cell
VlnPlot(immune.combined, features = 'Iqgap3',pt.size = 0)#IA

VlnPlot(immune.combined,features = c('Ptprc','Pecam1','Fibin','Igfbp7'),pt.size = 0)

#VlnPlot(immune.combined,features = c('Alb','Afp','Albumin','Acly','Apoa','Asl','Ass'),pt.size = 0) #hepatocyte
VlnPlot(immune.combined,features = c('Acta2','Col1a1','Hba-a1','Hba-a2','Hbb-bs'),pt.size = 0)#non-EC
VlnPlot(immune.combined,features = c('Cdh5','Pecam1','Prox1','Pdpn','Lyve-1'),pt.size = 0)#EC

VlnPlot(immune.combined,features = c('Clu','Plac8','Lrg1','Cd34'),pt.size = 0)#large artery
VlnPlot(immune.combined,features = c('Ly6c1','Rgcc','Timp4','Rbp7'),pt.size = 0)#large artery
VlnPlot(immune.combined,features = c('Efnb1','Glu1','Gja4','Sox17'),pt.size = 0)#capi arterial
VlnPlot(immune.combined,features = c('Atp13a3','Fbn1','Jag1','Adgrg6 ','Slc45a4'),pt.size = 0)#capi arterial
VlnPlot(immune.combined,features = c('Stab2','Kdr'),pt.size = 0)#capillary
VlnPlot(immune.combined,features = c('Mt1','Fcgr2b','Sepp1','Gatm'),pt.size = 0)#capillary
VlnPlot(immune.combined,features = c('Clec4g','Thbd'),pt.size = 0)#capi venous
VlnPlot(immune.combined,features = c('Rspo3','Slc26a1','Sema3d'),pt.size = 0)#capi venous
VlnPlot(immune.combined,features = c('Vwf','Edn1','Bgn'),pt.size = 0)#vein
VlnPlot(immune.combined,features = c('Alcam','Rbms3','Serpinb9b'),pt.size = 0)#vein
VlnPlot(immune.combined,features = c('Stmn1','Cdkn1a','Mki67','Top2a'),pt.size = 0)#proliferating
VlnPlot(immune.combined,features = c('Ccl21a','Prss23','Prox1'),pt.size = 0)#lymphatic
VlnPlot(immune.combined,features = c('Cxcl12','Efnb2','Rbp7'),pt.size = 0)#large artery
VlnPlot(immune.combined,features = c('Nr2f2','Nrp2','Ackr3','Vcam1'),pt.size = 0)
VlnPlot(immune.combined,features = c('Ctss','Cd74'),pt.size = 0)#large artery
VlnPlot(immune.combined,features = c('Vwa1','Aplnr'),pt.size = 0)
VlnPlot(immune.combined,features = 'Cd36',pt.size = 0)


cluster15.markers <- FindConservedMarkers(immune.combined,ident.1 = 15,grouping.var = 'stim')
VlnPlot(immune.combined,features = c('Cst3','Fabp4','Gm12840','Apold1','Cfh','Hspb1'),pt.size = 0)


immune.combined <- RenameIdents(immune.combined,'8'='artery','11'='artery','14'='artery',
                                '6'='cap.art','15'='cap.art',
                                '0'='capillary','1'='capillary','2'='capillary','3'='capillary','4'='capillary','16'='capillary','17'='capillary',
                                '7'='cap.ven',
                                '9'='vein',
                                '13'='lymphatic',
                                '5'='proliferating','10'='proliferating','12'='proliferating')
                                 
                                #'11'='IA cell',
                                #'13'='tip cell')

Idents(immune.combined) = immune.combined$seurat_clusters
new.cluster.ids <- c("capillary","capillary","capillary","capillary","capillary","capillary",'cap.art','cap.ven','artery',
                     'vein','proliferating','artery','proliferating','lymphatic','artery','cap.art','capillary','artery')
names(new.cluster.ids) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined, new.cluster.ids)
immune.combined$celltype = as.character(Idents(immune.combined))


DimPlot(immune.combined, label = TRUE)



cmarker = FindConservedMarkers(immune.combined,ident.1 = 'capillary',grouping.var = 'stim')
cvmarker = FindConservedMarkers(immune.combined,ident.1 = 'cap.ven',grouping.var = 'stim')
vmarker = FindConservedMarkers(immune.combined,ident.1 = 'vein',grouping.var = 'stim')
camarker = FindConservedMarkers(immune.combined,ident.1 = 'cap.art',grouping.var = 'stim')
lmarker = FindConservedMarkers(immune.combined,ident.1 = 'lymphatic',grouping.var = 'stim')
amarker = FindConservedMarkers(immune.combined,ident.1 = 'artery',grouping.var = 'stim')
tipmarker = FindConservedMarkers(immune.combined,ident.1 = 'tip cell',grouping.var = 'stim')
iamarker = FindConservedMarkers(immune.combined,ident.1 = 'IA',grouping.var = 'stim')

Capillary=subset(immune.combined,idents = 'capillary')
Cap.ven=subset(immune.combined,idents = 'cap.ven')
Vein=subset(immune.combined,idents = 'vein')
Cap.art=subset(immune.combined,idents = 'cap.art')
Proliferating=subset(immune.combined,idents = 'proliferating')
Lymphatic=subset(immune.combined,idents = 'lymphatic')
Artery=subset(immune.combined,idents = 'artery')

table(Idents(Artery))

table(Cap.ven$stim)  

prop.table(table(Idents(immune.combined), immune.combined$stim))

cell.prop<-as.data.frame(prop.table(table(Idents(immune.combined), immune.combined$stim)))
colnames(cell.prop)<-c("cluster","stim","proportion")
ggplot(cell.prop,aes(stim,proportion,fill=cluster))+
  geom_bar(stat="identity",position= "fill")+
  coord_flip()+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))

FeaturePlot(immune.combined, features = "Mki67", split.by = "stim", max.cutoff = 3, 
            cols = c("grey", "blue"))
FeaturePlot(immune.combined, features = "Vwf", split.by = "stim", max.cutoff = 3, 
            cols = c("grey", "red"))

FeaturePlot(immune.combined, features = c('Stmn1','Cdkn1a','Mki67','Top2a'))
FeaturePlot(immune.combined, features = 'Vwf')#vein
FeaturePlot(immune.combined, features = 'Vcam1')#capi venous
FeaturePlot(immune.combined, features = 'Ccl21a')#lymphatic
FeaturePlot(immune.combined, features = 'Efnb1')#capi arterial
FeaturePlot(immune.combined, features = 'Stab2')#capillary
FeaturePlot(immune.combined, features = 'Top2a')#proliferating
FeaturePlot(immune.combined, features = 'Cxcl12')#large artery

RidgePlot(immune.combined, features = "Mki67")

VlnPlot(immune.combined,features = 'Vwf',pt.size = 0)#vein
VlnPlot(immune.combined,features = 'Thbd',pt.size = 0)#cap ven
VlnPlot(immune.combined,features = 'Efnb1',pt.size = 0)#cap art
VlnPlot(immune.combined,features = 'Stab2',pt.size = 0)#capillary
VlnPlot(immune.combined,features = 'Cxcl12',pt.size = 0)#artery
VlnPlot(immune.combined,features = 'Top2a',pt.size = 0)#proliferating
VlnPlot(immune.combined,features = 'Ccl21a',pt.size = 0)#lymphatic

immunecombinedmarkers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
immunecombinedmarkers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top40 <- immunecombinedmarkers %>% group_by(cluster) %>% top_n(n = 40, wt = avg_logFC)
DoHeatmap(immune.combined, features = top50$gene) 




table(immune.combined$stim)

markers <- read.csv("zlycellinfo.csv")


?aggregate
aggdata <- aggregate(exprs$Mki67,list(exprs$celltype),mean)

?FetchData
#接findconservedmarkers后
cagene <- as.matrix(rownames(camarker))
cvgene <- as.matrix(rownames(cvmarker))
agene <- as.matrix(rownames(amarker))
cgene <- as.matrix(rownames(cmarker))
vgene <- as.matrix(rownames(vmarker))
pgene <- as.matrix(rownames(pmarker))
lgene <- as.matrix(rownames(lmarker))
iagene <- as.matrix(rownames(iamarker))
tipgene <- as.matrix(rownames(tipmarker))

exprs <- FetchData(immune.combined, vars  = as.character(c(cagene[1:20],cvgene[1:20],agene[1:20],cgene[1:20],
                                                      vgene[1:20],pgene[1:20],lgene[1:20],iagene[1:20],tipgene[1:20])[!duplicated(c(cagene[1:20],cvgene[1:20],agene[1:20],cgene[1:20],
                                                                                                         vgene[1:20],pgene[1:20],lgene[1:20],iagene[1:20],tipgene[1:20]))]))
#agene[1:20][!(agene[1:20]%in%cgene[1:20])]
celltype <- immune.combined@meta.data$celltype
exprs <- cbind(exprs,celltype)
aggdata <- aggregate(exprs[,1:147],list(exprs$celltype),mean)
rownames(aggdata) <- aggdata[,1]
aggdata <- subset.data.frame(aggdata[,2:148])
#注意这里调整行顺序
b <- subset.data.frame(aggdata[2:3,])
c <- subset.data.frame(aggdata[7,])
b <- rbind(b,c)
aggdata <- b
tagg <- t(aggdata)
tagg <- data.frame(tagg)

pheatmap(aggdata,scale = 'column',
         cluster_rows = 0,cluster_cols = 0,
         border_color = 0,
         show_colnames = T,
         show_rownames = T,
         color = colorRampPalette(colors = c("blue","white","red"))(100)
         )






#相关性分析
AverageExp<-AverageExpression(immune.combined,features=unique(top10$gene))
coorda<-corr.test(AverageExp$RNA,AverageExp$RNA,method="spearman")
pheatmap(coorda$r)

write.table(a,'adaptedagg.csv',row.names = T,col.names = T)
write.table(aggdata,'aggdata.csv',row.names = T,col.names = T)
write.table(immunecombinedmarkers,'immunecombinedmarkers.csv',row.names = T,col.names = T)
write.table(top10,'top10markers.csv',row.names = T,col.names = T)


devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')


RA_matrix<-as(as.matrix(immune.combined@assays$RNA@data), 'sparseMatrix')
feature_ann<-data.frame(gene_id=rownames(RA_matrix),gene_short_name=rownames(RA_matrix))
rownames(feature_ann)<-rownames(RA_matrix)
RA_fd<-new("AnnotatedDataFrame", data = feature_ann)
sample_ann<-immune.combined@meta.data
rownames(sample_ann)<-colnames(RA_matrix)
RA_pd<-new("AnnotatedDataFrame", data =sample_ann)

RA.cds<-newCellDataSet(RA_matrix,phenoData =RA_pd,featureData =RA_fd,expressionFamily=negbinomial.size())

BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
Heatmap(a,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = TRUE)
mark_gene <- c('Ly6a','Thbd','Ly6c1','Stab2','Rspo3','Mki67','Ccl21a')
gene_pos <- which(colnames(a) %in% mark_gene)
col_anno <-  columnAnnotation(mark_gene = anno_mark(at = gene_pos, 
                                                 labels = mark_gene))
Heatmap(a,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = TRUE,
        top_annotation = col_anno,
        column_title = NULL)

uncertain <- subset(immune.combined,idents = c(7,9))
DimPlot(uncertain,label = 1)
DefaultAssay(uncertain) <- "integrated"
uncertain <- ScaleData(uncertain, verbose = 1)
uncertain <- RunPCA(uncertain, npcs = 30, verbose = 1)
uncertain <- RunUMAP(uncertain, reduction = "pca", dims = 1:20,min.dist = 0.4)
uncertain <- FindNeighbors(uncertain, reduction = "pca", dims = 1:20)
uncertain <- FindClusters(uncertain, resolution = 0.7)
DimPlot(uncertain,label = 1)
DefaultAssay(uncertain) <- "RNA"
VlnPlot(uncertain,features = c('Bmp4','Prss23','Tgfb2','Cyp1b1','Alcam','Vwf','Rspo3','Entpd1','Ctsh','Wnt9b'),pt.size = 0)#vein
VlnPlot(uncertain,features = c(#'Lhx6','Selp','Wnt2','Thbd','Cdh13','Cd9','Slc26a10','Jam2','Cldn5','Bcam',
                               'App','Cd36',"Flrt2","Lbp","Ptgs1","Ackr4","Plvap",'Slco3a1','Lrg1','Ramp3','Qsox1'),pt.size = 0)#cap.ven
cluster1.markers <- FindMarkers(uncertain, ident.1 = 1, min.pct = 0)
cluster0.markers <- FindMarkers(uncertain, ident.1 = 0, min.pct = 0.25)
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
VlnPlot(uncertain,features = c('Clec4g','Wdr89','Smagp','S100a1','Cldn5','Gpihbp1','Rps12-ps3','Gm11808','Dnase1l3','Plpp1'),pt.size = 0)
uncertain <- RenameIdents(uncertain,'0'='capven','3'='capven','4'='capven','6'='capven',
                          '1'='vein','2'='vein','5'='vein')
DimPlot(uncertain,label = 1)
uncertain$celltype = as.character(Idents(uncertain))
cvmarker = FindConservedMarkers(uncertain,ident.1 = 'capven',grouping.var = 'stim')
vmarker = FindConservedMarkers(uncertain,ident.1 = 'vein',grouping.var = 'stim')
cvgene <- as.matrix(rownames(cvmarker))
vgene <- as.matrix(rownames(vmarker))

exprs <- FetchData(uncertain, vars  = as.character(c(cvgene,vgene)))
celltype <- uncertain@meta.data$celltype
exprs <- cbind(exprs,celltype)
aggdata <- aggregate(exprs[,1:289],list(exprs$celltype),mean)
rownames(aggdata) <- aggdata[,1]
aggdata <- subset.data.frame(aggdata[,2:289])
FeaturePlot(uncertain, features = 'Vwf')#Ctsh,Bmp4,Edn1
FeaturePlot(uncertain, features = 'Dnase1l3')#Gatm,Plxnc1,Maf,Mrc1

tagg <- t(aggdata)
tagg <- data.frame(tagg)
targets <- as.matrix(tagg)
group <- rep(c('cv', 'v'), each = 1)
head(targets)
dgelist <- DGEList(counts = targets, group = group)
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, ,keep.lib.sizes = FALSE]
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
#plotMDS(dgelist_norm, col = rep(c('red', 'blue'), each = 1), dim = c(1, 2))
design <- model.matrix(~group)  
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
plotBCV(dge)

FeaturePlot(immune.combined, features = 'Top2a',split.by = 'stim')
FeaturePlot(immune.combined, features = 'Foxo3')
VlnPlot(immune.combined, features = 'Top2a',pt.size = 0,split.by = 'stim')
VlnPlot(immune.combined, features = 'Ube2c',pt.size = 0)

DefaultAssay(Proliferating) <- "integrated"
Proliferating <- ScaleData(Proliferating, verbose = 1)
Proliferating <- RunPCA(Proliferating, npcs = 30, verbose = 1)
Proliferating <- RunUMAP(Proliferating, reduction = "pca", dims = 1:20,min.dist = 0.4)
Proliferating <- FindNeighbors(Proliferating, reduction = "pca", dims = 1:20)
Proliferating <- FindClusters(Proliferating, resolution = 0.5)
DimPlot(Proliferating,label = 1)
DefaultAssay(Proliferating) <- "RNA"
VlnPlot(Proliferating, features = 'Adm',pt.size = 0)#tip cell
VlnPlot(Proliferating, features = 'Iqgap3',pt.size = 0)#IA
VlnPlot(Proliferating, features = 'Cenpe',pt.size = 0)#tip
VlnPlot(Proliferating, features = 'Mki67',pt.size = 0)
Proliferating <- RenameIdents(Proliferating,'0'='tipcell','4'='tipcell','1'='IA','3'='tipcell','2'='IA')
Proliferating$celltype = as.character(Idents(Proliferating))
DimPlot(Proliferating,label = TRUE)
table(subset(Proliferating,subset=celltype== 'IA')$stim)
table(subset(Proliferating,subset=celltype== 'tipcell')$stim)

tipmarker = FindConservedMarkers(Proliferating,ident.1 = 'tip cell',grouping.var = 'stim')
iamarker = FindConservedMarkers(Proliferating,ident.1 = 'IA',grouping.var = 'stim')
FeaturePlot(Proliferating, features = 'Top2a')
VlnPlot(immune.combined, features = 'Adm',pt.size = 0)

Proliferating = subset(Proliferating,subset = stim!='wt-ph0d')

FeaturePlot(Proliferating, features = c('Adm','Cenpe','Iqgap3'),pt.size = 0)
VlnPlot(Proliferating, features = 'Gas6',pt.size = 0.1,split.by = 'stim')
tva <- FindMarkers(Proliferating, ident.1 = 'tipcell', ident.2 = 'IA', min.pct = 0.2,logfc.threshold = 0.2)


plots <- VlnPlot(Proliferating, features = 'Gas6', split.by = "stim", group.by = "celltype", 
                 pt.size = 0.1, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

prop.table(table(Idents(Proliferating), Proliferating$stim))
table(subset(Proliferating,subset=celltype== 'IA')$stim)/sum(table(Proliferating$stim))
