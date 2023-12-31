# Sergio Alías, 20231027
# Last modified 20231027

# Script for integrate neuron organoid single-cell data (two genotype conditions)

library(Seurat)
library(DoubletFinder)
library(scCustomize)
library(patchwork)
library(harmony)

dir.create("integration_outs")

dataset_loc <- "/mnt2/fscratch/users/bio_267_uma/sergioalias/NGS_projects/NEURORG/count_results"
ids.WT <- c("SAMPLE_B",
            "SAMPLE_C",
            "SAMPLE_E",
            "SAMPLE_H",
            "SAMPLE_I",
            "SAMPLE_K")

ids.mutant <- c("SAMPLE_A",
                "SAMPLE_D",
                "SAMPLE_F",
                "SAMPLE_G",
                "SAMPLE_J",
                "SAMPLE_L")

seu.list.WT <- sapply(ids.WT, function(i){ # Loading
  d10x <- Read10X(file.path(dataset_loc,i,"cellranger_0000",i,"outs/filtered_feature_bc_matrix"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  seu <- CreateSeuratObject(counts = d10x, project = "neurorg", min.cells = 1, min.features = 1)
  seu@meta.data$samplename <- c(rep(i, nrow(seu@meta.data)))
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu <- subset(seu, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 20)
  seu
})

seu.list.mutant <- sapply(ids.mutant, function(i){ # Loading
  d10x <- Read10X(file.path(dataset_loc,i,"cellranger_0000",i,"outs/filtered_feature_bc_matrix"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  seu <- CreateSeuratObject(counts = d10x, project = "neurorg", min.cells = 1, min.features = 1)
  seu@meta.data$samplename <- c(rep(i, nrow(seu@meta.data)))
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu <- subset(seu, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 20)
  seu
})


seu.list.WT$SAMPLE_B@meta.data$time <- c(rep("48", nrow(seu.list.WT$SAMPLE_B@meta.data)))
seu.list.WT$SAMPLE_C@meta.data$time <- c(rep("48", nrow(seu.list.WT$SAMPLE_C@meta.data)))
seu.list.WT$SAMPLE_E@meta.data$time <- c(rep("48", nrow(seu.list.WT$SAMPLE_E@meta.data)))
seu.list.WT$SAMPLE_H@meta.data$time <- c(rep("62", nrow(seu.list.WT$SAMPLE_H@meta.data)))
seu.list.WT$SAMPLE_I@meta.data$time <- c(rep("62", nrow(seu.list.WT$SAMPLE_I@meta.data)))
seu.list.WT$SAMPLE_K@meta.data$time <- c(rep("62", nrow(seu.list.WT$SAMPLE_K@meta.data)))


seu.list.mutant$SAMPLE_A@meta.data$time <- c(rep("48", nrow(seu.list.mutant$SAMPLE_A@meta.data)))
seu.list.mutant$SAMPLE_D@meta.data$time <- c(rep("48", nrow(seu.list.mutant$SAMPLE_D@meta.data)))
seu.list.mutant$SAMPLE_F@meta.data$time <- c(rep("48", nrow(seu.list.mutant$SAMPLE_F@meta.data)))
seu.list.mutant$SAMPLE_G@meta.data$time <- c(rep("62", nrow(seu.list.mutant$SAMPLE_G@meta.data)))
seu.list.mutant$SAMPLE_J@meta.data$time <- c(rep("62", nrow(seu.list.mutant$SAMPLE_J@meta.data)))
seu.list.mutant$SAMPLE_L@meta.data$time <- c(rep("62", nrow(seu.list.mutant$SAMPLE_L@meta.data)))


seu.WT <- Merge_Seurat_List(list_seurat = seu.list.WT)
seu.WT <- NormalizeData(seu.WT)
seu.WT <- FindVariableFeatures(seu.WT, selection.method = "vst", nfeatures = 2000)
seu.WT <- ScaleData(seu.WT)
seu.WT <- RunPCA(seu.WT)

seu.WT <- RunHarmony(seu.WT, "samplename", plot_convergence = FALSE)
seu.WT <- RunUMAP(seu.WT, dims = 1:10, reduction = "harmony")
seu.WT <- FindNeighbors(seu.WT, dims = 1:10, reduction = "harmony")
seu.WT <- FindClusters(seu.WT, resolution = 0.6)

seu.markers.WT <- FindAllMarkers(seu.WT, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


seu.mutant <- Merge_Seurat_List(list_seurat = seu.list.mutant)
seu.mutant <- NormalizeData(seu.mutant)
seu.mutant <- FindVariableFeatures(seu.mutant, selection.method = "vst", nfeatures = 2000)
seu.mutant <- ScaleData(seu.mutant)
seu.mutant <- RunPCA(seu.mutant)

seu.mutant <- RunHarmony(seu.mutant, "samplename", plot_convergence = FALSE)
seu.mutant <- RunUMAP(seu.mutant, dims = 1:10, reduction = "harmony")
seu.mutant <- FindNeighbors(seu.mutant, dims = 1:10, reduction = "harmony")
seu.mutant <- FindClusters(seu.mutant, resolution = 0.6)

seu.markers.mutant <- FindAllMarkers(seu.mutant, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Figures

setwd("outs")

saveRDS(seu.WT, file = "seu.WT.RDS")
saveRDS(seu.mutant, file = "seu.mutant.RDS")
saveRDS(seu.markers.WT, file = "seu.markers.WT.RDS")
saveRDS(seu.markers.mutant, file = "seu.markers.mutant.RDS")