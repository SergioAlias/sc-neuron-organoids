# Sergio Alías, 20231214
# Last modified 20231214

# Mini script for adding manual annotations to integrated Seurat objects

library(Seurat)

seu.ints <- c("48",
              "62",
              "WT",
              "mutant")

for (i in seu.ints){
  seu <- readRDS(file.path("outs", paste0("seu.", i, ".RDS")))
  new.cluster.ids <- read.table(file = paste0("cluster-celltypes-neurorg-", i),
                                sep = '\t',
                                header = FALSE)[,2]
  names(new.cluster.ids) <- levels(seu)
  seu <- RenameIdents(seu, new.cluster.ids)
  saveRDS(seu, file = file.path("outs", paste0("seu.", i, ".RDS")))
}