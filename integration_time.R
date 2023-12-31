# Sergio Alías, 20231025
# Last modified 20231025

# Script for integrate neuron organoid single-cell data (two time conditions)

library(Seurat)
library(DoubletFinder)
library(scCustomize)
library(patchwork)
library(harmony)

dir.create("integration_outs")

dataset_loc <- "/mnt2/fscratch/users/bio_267_uma/sergioalias/NGS_projects/NEURORG/count_results"
ids.48 <- c("SAMPLE_A",
            "SAMPLE_B",
            "SAMPLE_C",
            "SAMPLE_D",
            "SAMPLE_E",
            "SAMPLE_F")

ids.62 <- c("SAMPLE_G",
            "SAMPLE_H",
            "SAMPLE_I",
            "SAMPLE_J",
            "SAMPLE_K",
            "SAMPLE_L")

seu.list.48 <- sapply(ids.48, function(i){ # Loading
  d10x <- Read10X(file.path(dataset_loc,i,"cellranger_0000",i,"outs/filtered_feature_bc_matrix"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  seu <- CreateSeuratObject(counts = d10x, project = "neurorg", min.cells = 1, min.features = 1)
  seu@meta.data$samplename <- c(rep(i, nrow(seu@meta.data)))
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu <- subset(seu, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 20)
  seu
})

seu.list.62 <- sapply(ids.62, function(i){ # Loading
  d10x <- Read10X(file.path(dataset_loc,i,"cellranger_0000",i,"outs/filtered_feature_bc_matrix"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  seu <- CreateSeuratObject(counts = d10x, project = "neurorg", min.cells = 1, min.features = 1)
  seu@meta.data$samplename <- c(rep(i, nrow(seu@meta.data)))
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu <- subset(seu, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 20)
  seu
})


seu.list.48$SAMPLE_A@meta.data$genotype <- c(rep("mutant", nrow(seu.list.48$SAMPLE_A@meta.data)))
seu.list.48$SAMPLE_B@meta.data$genotype <- c(rep("WT", nrow(seu.list.48$SAMPLE_B@meta.data)))
seu.list.48$SAMPLE_C@meta.data$genotype <- c(rep("WT", nrow(seu.list.48$SAMPLE_C@meta.data)))
seu.list.48$SAMPLE_D@meta.data$genotype <- c(rep("mutant", nrow(seu.list.48$SAMPLE_D@meta.data)))
seu.list.48$SAMPLE_E@meta.data$genotype <- c(rep("WT", nrow(seu.list.48$SAMPLE_E@meta.data)))
seu.list.48$SAMPLE_F@meta.data$genotype <- c(rep("mutant", nrow(seu.list.48$SAMPLE_F@meta.data)))


seu.list.62$SAMPLE_G@meta.data$genotype <- c(rep("mutant", nrow(seu.list.62$SAMPLE_G@meta.data)))
seu.list.62$SAMPLE_H@meta.data$genotype <- c(rep("WT", nrow(seu.list.62$SAMPLE_H@meta.data)))
seu.list.62$SAMPLE_I@meta.data$genotype <- c(rep("WT", nrow(seu.list.62$SAMPLE_I@meta.data)))
seu.list.62$SAMPLE_J@meta.data$genotype <- c(rep("mutant", nrow(seu.list.62$SAMPLE_J@meta.data)))
seu.list.62$SAMPLE_K@meta.data$genotype <- c(rep("WT", nrow(seu.list.62$SAMPLE_K@meta.data)))
seu.list.62$SAMPLE_L@meta.data$genotype <- c(rep("mutant", nrow(seu.list.62$SAMPLE_L@meta.data)))


seu.48 <- Merge_Seurat_List(list_seurat = seu.list.48)
seu.48 <- NormalizeData(seu.48)
seu.48 <- FindVariableFeatures(seu.48, selection.method = "vst", nfeatures = 2000)
seu.48 <- ScaleData(seu.48)
seu.48 <- RunPCA(seu.48)

seu.48 <- RunHarmony(seu.48, "samplename", plot_convergence = FALSE)
seu.48 <- RunUMAP(seu.48, dims = 1:10, reduction = "harmony")
seu.48 <- FindNeighbors(seu.48, dims = 1:10, reduction = "harmony")
seu.48 <- FindClusters(seu.48, resolution = 0.6)

seu.markers.48 <- FindAllMarkers(seu.48, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


seu.62 <- Merge_Seurat_List(list_seurat = seu.list.62)
seu.62 <- NormalizeData(seu.62)
seu.62 <- FindVariableFeatures(seu.62, selection.method = "vst", nfeatures = 2000)
seu.62 <- ScaleData(seu.62)
seu.62 <- RunPCA(seu.62)

seu.62 <- RunHarmony(seu.62, "samplename", plot_convergence = FALSE)
seu.62 <- RunUMAP(seu.62, dims = 1:10, reduction = "harmony")
seu.62 <- FindNeighbors(seu.62, dims = 1:10, reduction = "harmony")
seu.62 <- FindClusters(seu.62, resolution = 0.6)

seu.markers.62 <- FindAllMarkers(seu.62, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Figures

setwd("outs")

saveRDS(seu.48, file = "seu.48.RDS")
saveRDS(seu.62, file = "seu.62.RDS")
saveRDS(seu.markers.48, file = "seu.markers.48.RDS")
saveRDS(seu.markers.62, file = "seu.markers.62.RDS")