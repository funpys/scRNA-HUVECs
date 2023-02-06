library(DropletUtils)
library(scater)
load("ensemblGenes2017-08-09.RData")


### load raw data
sce1 <- read10xCounts("./data_3d/w-LF-GEX/")
colnames(sce1) = paste0("w-",sce1$Barcode)

sce2 <- read10xCounts("./data_3d/wo-LF-GEX/")
colnames(sce2) = paste0("wo-",sce2$Barcode)

sce = cbind(sce1, sce2)

### rank plot
minColSum=160
countGBM = counts(sce)
countGBMFiltered = countGBM[,colSums(countGBM)>160]
dim(countGBMFiltered)

count = colSums(countGBMFiltered)
pdf("HUVEC_3D_rankplot.pdf")
plot(1:length(count), sort(count, decreasing=T), pch=20, log="y")
abline(h=1000, col="red", lty=2, lwd=2)
dev.off()





### Create SingleCellExperiment object
sce = SingleCellExperiment(assays=list(counts=countGBMFiltered))
mtGenes = ensemblGenes$ensembl_gene_id[grepl("^MT-",ensemblGenes$external_gene_name)]
is.mito = rownames(sce) %in% mtGenes
length(is.mito[is.mito == TRUE])
sce_drop_filtered <- addPerCellQC(sce, subsets=list(Mito=is.mito))
sce_drop_filtered = runColDataPCA(sce_drop_filtered, variables=list(
  "sum", "detected", "subsets_Mito_percent", "subsets_Mito_sum", "subsets_Mito_detected", "subsets_Mito_percent"
))
sce_drop_filtered$log10_total_counts = log10(sce_drop_filtered$sum)

### total_counts
log_count_cutoff = 3

png(paste0("HUVECs_3d_wo_histogram_totalcount.png"))
h <- hist(
  sce_drop_filtered$log10_total_counts,
  breaks = 100
)
cuts <- cut(h$breaks, c(0,log_count_cutoff,Inf))
plot(h, col=c("red","white")[cuts])
abline(v = log_count_cutoff, col="red")
dev.off()

png(paste0("HUVECs_3d_wo_PCA_totalcount.png"))
p=plotReducedDim(sce_drop_filtered,
                 dimred = "PCA_coldata",
                 colour_by = "log10_total_counts"
)
print(p)
dev.off()

sce_drop_filtered$filter_counts = sce_drop_filtered$log10_total_counts<log_count_cutoff
# sum(sce_drop_filtered$log10_total_counts>log_count_cutoff)

png(paste0("HUVECs_3d_wo_PCA_totalcount_filter.png"))
p=plotReducedDim(sce_drop_filtered,
                 dimred = "PCA_coldata",
                 colour_by = "filter_counts"
)
print(p)
dev.off()



### pct_counts_Mt
Mt_ratio_cutoff = 10

png(paste0("HUVECs_3d_wo_histogram_MT_ratio.png"))
h <- hist(
  sce_drop_filtered$subsets_Mito_percent,
  breaks = 100
)
cuts <- cut(h$breaks, c(0,Mt_ratio_cutoff,Inf))
plot(h, col=c("white", "red")[cuts])
abline(v = Mt_ratio_cutoff, col="red")
dev.off()

# sum(sce_drop_filtered$pct_counts_Mt<Mt_ratio_cutoff)
# sce_drop_filtered$pct_counts_Mt[sce_drop_filtered$pct_counts_Mt>10] = 10

png(paste0("HUVECs_3d_wo_PCA_MT_ratio.png"))
p=plotReducedDim(sce_drop_filtered,
                 dimred = "PCA_coldata",
                 colour_by = "subsets_Mito_percent"
)
print(p)
dev.off()

sce_drop_filtered$filter_MT = sce_drop_filtered$subsets_Mito_percent>Mt_ratio_cutoff
sum(sce_drop_filtered$subsets_Mito_percent<Mt_ratio_cutoff)

png(paste0("HUVECs_3d_wo_PCA_MT_ratio_filter.png"))
p=plotReducedDim(sce_drop_filtered,
                 dimred = "PCA_coldata",
                 colour_by = "filter_MT"
)
print(p)
dev.off()



### filter low-quality cells
sce_drop_filtered$use <- (
  !sce_drop_filtered$filter_counts &
    !sce_drop_filtered$filter_MT 
)
sce_drop_filtered$use = !sce_drop_filtered$use
sum(!sce_drop_filtered$use)

png(paste0("HUVECs_3d_wo_PCA_filter.png"))
p=plotReducedDim(sce_drop_filtered,
                 dimred = "PCA_coldata",
                 colour_by = "use"
)
print(p)
dev.off()

sce_QC = sce_drop_filtered[ , !sce_drop_filtered$use]
save(sce_QC, file = paste0("sce_drop_filtered_QC_HUVECs_3d.RData"))
