library(scran)
library(Seurat)
library(dplyr)
library(magrittr)
load("sce_drop_filtered_QC_HUVECs_3d.RData")

### Normalization
clusters <- quickCluster(sce_QC)
sce_QC <- computeSumFactors(sce_QC, cluster=clusters)
summary(sizeFactors(sce_QC))
sce_norm = logNormCounts(sce_QC)


### Create Seurat object
HUVEC_3d_seurat = CreateSeuratObject(counts = counts(sce_norm), project = "HUVEC_3d")
HUVEC_3d_seurat@assays$RNA@data = logcounts(sce_norm)
HUVEC_3d_seurat$MT_ratio = sce_drop_filtered$subsets_Mito_percent

### Get HVGs
dec <- modelGeneVar(sce_norm)
table(dec$FDR <0.05)[2]
top.hvgs <- getTopHVGs(dec, fdr.threshold = 0.05)
table(dec$FDR <0.05)
HUVEC_3d_seurat@assays$RNA@var.features = top.hvgs

### Get cell cycle score & correct cell cycle effect
cc.genes <- readLines(con = "D:/R_code/PC9_code/regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
S.genes <- ensemblGenes[1][match(s.genes, ensemblGenes[,2],nomatch = 0),]
G2M.genes <- ensemblGenes[1][match(g2m.genes, ensemblGenes[,2],nomatch = 0),]
HUVEC_3d_seurat <- CellCycleScoring(HUVEC_3d_seurat, s.features = S.genes, g2m.features = G2M.genes, set.ident = F)
HUVEC_3d_seurat = ScaleData(HUVEC_3d_seurat, vars.to.regress = c("S.Score", "G2M.Score"))

### Visualization
HUVEC_3d_seurat <- RunPCA(HUVEC_3d_seurat)
ElbowPlot(HUVEC_3d_seurat, ndims = 50)
dim = 1:30
HUVEC_3d_seurat <- RunUMAP(HUVEC_3d_seurat, dims = dim)
HUVEC_3d_seurat <- FindNeighbors(HUVEC_3d_seurat, dims = dim)
HUVEC_3d_seurat <- FindClusters(HUVEC_3d_seurat)
HUVEC_3d_seurat$condition = strsplit(colnames(HUVEC_3d_seurat), "_") %>% sapply(extract2, 1)

DimPlot(HUVEC_3d_seurat, reduction = "umap", label = T, label.size = 5, pt.size = .8)
DimPlot(HUVEC_3d_seurat, reduction = "umap", label = F, label.size = 5, pt.size = .8, group.by = "condition", cols = c("#e41a1c", "#377eb8"))
DimPlot(HUVEC_3d_seurat, reduction = "umap", label = T, label.size = 5, pt.size = .8, group.by = "Phase")
DimPlot(HUVEC_3d_seurat, reduction = "umap", label = T, label.size = 5, pt.size = .8, group.by = "cell_type")



### Filter low quality cells
HUVEC_3d_seurat$nCount_RNA_log10 = log10(HUVEC_3d_seurat$nCount_RNA)
VlnPlot(HUVEC_3d_seurat, pt.size = 0, features = "nCount_RNA_log10")
VlnPlot(HUVEC_3d_seurat, pt.size = 0, features = "MT_ratio")

HUVEC_3d_seurat = subset(HUVEC_3d_seurat, cells = unique(c(cells,colnames(HUVEC_3d_seurat)[HUVEC_3d_seurat@active.ident %in% c(2,9,10,4,
                                                                                                                               6,5,11,14,15,3,7,1,13)])))
### Normalization
sce = SingleCellExperiment(assays=list(counts=HUVEC_3d_seurat@assays$RNA@counts))
sce_QC = sce[rowSums(counts(sce))>0,]
clusters <- quickCluster(sce_QC)
sce_QC <- computeSumFactors(sce_QC, cluster=clusters)
summary(sizeFactors(sce_QC))
sce_norm = logNormCounts(sce_QC)
HUVEC_3d_seurat@assays$RNA@data = logcounts(sce_norm)

### Get HVGs
dec <- modelGeneVar(sce_norm)
table(dec$FDR <0.05)[2]
# top.hvgs <- getTopHVGs(dec, n=2000)
top.hvgs <- getTopHVGs(dec, fdr.threshold = 0.05)
table(dec$FDR <0.05)
HUVEC_3d_seurat@assays$RNA@var.features = top.hvgs

### Visualization
HUVEC_3d_seurat = ScaleData(HUVEC_3d_seurat, vars.to.regress = c("S.Score", "G2M.Score"))
HUVEC_3d_seurat <- RunPCA(HUVEC_3d_seurat)

ElbowPlot(HUVEC_3d_seurat, ndims = 50)

dim = 1:30
HUVEC_3d_seurat <- RunUMAP(HUVEC_3d_seurat, dims = dim)
HUVEC_3d_seurat <- FindNeighbors(HUVEC_3d_seurat, dims = dim)
HUVEC_3d_seurat <- FindClusters(HUVEC_3d_seurat)

DimPlot(HUVEC_3d_seurat, reduction = "umap", label = F, label.size = 5, pt.size = .8, group.by = "condition", cols = c("#e41a1c", "#377eb8"))
DimPlot(HUVEC_3d_seurat, reduction = "umap", label = T, label.size = 5, pt.size = .8, group.by = "Phase")
DimPlot(HUVEC_3d_seurat, reduction = "umap", label = T, label.size = 5, pt.size = .8, group.by = "cell_type")
DimPlot(HUVEC_3d_seurat, reduction = "umap", label = T, label.size = 5, pt.size = .8, group.by = "filtered")
