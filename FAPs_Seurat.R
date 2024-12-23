library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(sctransform)
library(harmony)
library(Matrix)
library(future)
library(pheatmap)

GGPP.data <- Read10X(data.dir = "/~")
Statin.data <- Read10X(data.dir = "/~")
WT.data <- Read10X(data.dir = "/~")

GGPP <- CreateSeuratObject(counts = GGPP.data, project = "GGPP_co", min.cells = 3, min.features = 200)
Statin <- CreateSeuratObject(counts = Statin.data, project = "Statin", min.cells = 3, min.features = 200)
WT <- CreateSeuratObject(counts = WT.data, project = "WT_FAP", min.cells = 3, min.features = 200)

GGPP[["percent.mt"]] <- PercentageFeatureSet(GGPP, pattern = "^MT-")
Statin[["percent.mt"]] <- PercentageFeatureSet(Statin, pattern = "^MT-")
WT[["percent.mt"]] <- PercentageFeatureSet(WT, pattern = "^MT-")

VlnPlot(GGPP, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Statin, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(WT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

GGPP <- subset(GGPP, subset = percent.mt < 3)
Statin <- subset(Statin, subset = percent.mt < 15)
WT <- subset(WT, subset = percent.mt < 5)

combined <- merge(WT, y = c(Statin, GGPP), add.cell.ids = c("WT_FAP", "Statin", "GGPP"), project = "FAPS")

DefaultAssay(combined) <- "RNA"
combined1 <- SCTransform(combined, vars.to.regress = "percent.mt")
combined2 <- RunPCA(combined1, npcs = 50)
combined2 <- RunHarmony(combined2, group.by.vars = "orig.ident", dims = 1:10)
combined2 <- FindNeighbors(combined2, reduction = "harmony", dims = 1:10)
combined2 <- FindClusters(combined2, resolution = 0.2)
combined2 <- RunUMAP(combined2, reduction = "harmony", dims = 1:10)

DimPlot(combined2, reduction = "umap", group.by = c("seurat_clusters"), pt.size = 0.3,label = F, label.size = 4)

combined2 <- PrepSCTFindMarkers(combined2)
combined.markers <- FindAllMarkers(combined2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

combined2$seurat_clusters <- factor(combined2$seurat_clusters, levels = c("FAPs 1","FAPs 2","FAPs 3","FAPs 4","FAPs 5","CD90+ FAPs"))