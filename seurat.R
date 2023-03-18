# Following this tut: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html


library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = r"(C:\Users\Oscar Wright\Documents\youtube\RShiny_dashboard\pbmc3k_filtered_gene_bc_matrices\filtered_gene_bc_matrices\hg19)")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

pbmc <- RunUMAP(pbmc, dims = 1:10)

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
    "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc$celltype <- Idents(pbmc)


# Randomly assign cells: Age, Sex, sampleID
set.seed(123) # for reproducibility
sex <- sample(c("male", "female"), ncol(pbmc), replace = TRUE)
age <- round(runif(ncol(pbmc), 18, 65))
sampleID <- paste0("sample", sample(1:10, ncol(pbmc), replace = TRUE))

pbmc$sex <- sex
pbmc$age <- age
pbmc$sampleID <- sampleID

saveRDS(pbmc, file = r"(C:\Users\Oscar Wright\Documents\youtube\RShiny_dashboard\seurat_object.rds)")
