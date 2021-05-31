# 1. mtx file
# Load libraries
library(Seurat)
library(tidyverse)
library(SingleCellExperiment)
library(Matrix)
library(AnnotationHub)
library(ensembldb)
library(cowplot)
library(ggplot2)
library(scales)

# Read in `matrix.mtx`
x=readMM("~/matrix/raw_matrix.mtx.gz")

# Read in `genes.tsv`
gene <- read_tsv("GSE127465_gene_names_human_41861.tsv.gz", col_names = FALSE)
gene_ids <- gene$X1

# Read in `barcodes.tsv`
meta <- read_tsv("GSE127465_human_cell_metadata_54773x25.tsv.gz", col_names = F)
colnames(meta) <- meta[1,]
meta <- meta[-1,]
cell_ids <- meta$Barcode

# Create a sparse matrix for more efficient computation
x <- as(x, "dgCMatrix")
x[1:4,1:4]
dim(x)
x <- t(x)
x[1:4,1:4]
dim(x)

# Rename, row should be genes and column should be cell ids
colnames(x) <- cell_ids
rownames(x) <- gene_ids

# Creat meta data frame
meta=as.data.frame(cell_ids)
colnames(meta)=c('cell.name')
rownames(meta)=colnames(x)
# > rownames(meta)=colnames(x)
# Error in `.rowNamesDF<-`(x, value = value) : 
#   duplicate 'row.names' are not allowed
# In addition: Warning message:
#   non-unique values when setting 'row.names': 'bcAAAO', 'bcAAES', 'bcAAEX', 'bcAAFW', 'bcAAGF', 'bcAAHF', 'bcAAIH', 'bcAAKA', 'bcAAKU', 'bcAALE', 'bcAAMG', 'bcAAMK', 'bcAAMS', 'bcAANG', 'bcAAOQ', 'bcAAPP', 'bcAARI', 'bcAASS', 'bcAAUY', 'bcAAVQ', 'bcAAVX', 'bcAAWG', 'bcAAWH', 'bcAAWP', 'bcAAWX', 'bcAAZM', 'bcABAG', 'bcABBB', 'bcABBJ', 'bcABBM', 'bcABBR', 'bcABCF', 'bcABCM', 'bcABDA', 'bcABFP', 'bcABFR', 'bcABGR', 'bcABGW', 'bcABHL', 'bcABIS', 'bcABLX', 'bcABMK', 'bcABNK', 'bcABNT', 'bcABOK', 'bcABOX', 'bcABOZ', 'bcABQF', 'bcABQZ', 'bcABRS', 'bcABRT', 'bcABSW', 'bcABSX', 'bcABXQ', 'bcABZZ', 'bcACAO', 'bcACBE', 'bcACBG', 'bcACBJ', 'bcACCN', 'bcACDG', 'bcACDZ', 'bcACER', 'bcACFC', 'bcACFJ', 'bcACFL', 'bcACGA', 'bcACGE', 'bcACGF', 'bcACHD', 'bcACHS', 'bcACIU', 'bcACJQ', 'bcACKC', 'bcACKO', 'bcACMM', 'bcACND', 'bcACNL', 'bcACNM', 'bcACOB', 'bcACOD', 'bcACPR', 'bcACRM', 'bcACSI', 'bcACSZ', 'bcACTF', 'bcACUE', 'bcACUK', 'bcACVK', 'bcACWJ', 'bcACWM', 'bcACXO', 'bcACYL', 'bcACZI', 'bcADDV', 'bcADE [... truncated] 
colnames(x) <- as.character(1:54773)
rownames(meta) <- colnames(x)
head(meta)

# Creat Seurat object
object <- CreateSeuratObject(counts = x, 
                             meta=meta,
                             min.cells = 3, 
                             min.features = 200, 
                             project = "GSE127465")
object
object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# 2. csv file
library(data.table)
a=fread('~/huangzhen/GSM3148585_BC01_BLOOD1_counts.csv.gz',header = T)
a$V1 <- as.character(a$V1)
a<- t(a)
colnames(a)<- a[1,]
a <- a[-1,]
v1 <- rownames(a)
v1 <- as.data.frame(v1)
colnames(v1)<- "V1"
a <- cbind(v1,a)
a<-data.frame(a,row.names = NULL)
length(a$V1)
length(unique(a$V1))
a1 <- a[!duplicated(a$V1),]
hg=a1$V1 #genes
dat=a1[,2:ncol(a1)]
rownames(dat)=hg
hg[grepl('^MT-',hg)]
colnames(dat)
rownames(dat)
meta=as.data.frame(colnames(dat))
colnames(meta)=c('cell.name')
rownames(meta)=colnames(dat)
head(meta)
## 前面大量的代码，都是数据预处理
library(Seurat)
dat[1:4,1:4]
class(dat)
dat[is.na(dat)] <- 0
# 重点是构建 Seurat对象
obj  <- CreateSeuratObject(counts = dat,
                           meta.data = meta,
                           min.cells = 3, min.features = 200,project = '10x_PBMC')

obj
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# 3. to create 10X files
#data <- download.file("https://links.jianshu.com/go?to=https%3A%2F%2Fwww.ncbi.nlm.nih.gov%2Fgeo%2Fdownload%2F%3Facc%3DGSM3270887%26format%3Dfile%26file%3DGSM3270887%255FcountTable%255FcolonCreMin%252Etxt%252Egz",destfile = "matrix.txt.gz",quiet=FALSE,cacheOK=TRUE,extra=getOption("download.file.extra"))
library(Matrix)
colon.data <- read.csv(file='GSM3270887_countTable_colonCreMin.txt.gz', sep="\t", header = T, row.names = 1)
colon.data <- Matrix(as.matrix(colon.data), sparse=T)

ngenes <- nrow(colon.data)
psedu_gene.ids <- paste0("ENSG0000", seq_len(ngenes))

writeMM(obj = colon.data, file="./matrix.mtx")
write.table(data.frame(psedu_gene.ids,rownames(colon.data)), file="./genes.tsv", 
            col.names=F,row.names = F, sep = "\t", quote=FALSE)
write(x = colnames(colon.data), file = "./barcodes.tsv")

data <- Read10X(data.dir = "~/matrix/test/")
data <- CreateSeuratObject(counts = data, project = "colon", min.cells = 3, min.features = 200)

data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")#mus musculus
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
