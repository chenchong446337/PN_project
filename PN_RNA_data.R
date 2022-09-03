## for RNA data in PN

library(dplyr)
library(Seurat)
library(patchwork)

countsData<-read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/gene_data/PG_counts.csv", 
                     header = T, row.names=1)
pbmc <- CreateSeuratObject(counts = countsData, project = "Pn_single_cell", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

## normalize data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

## data scaling
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)


FeaturePlot(pbmc, features = c( "Oprk1", "Syt7"))

DimPlot(pbmc, reduction = "pca")

# umap dimension reduction-----
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

VlnPlot(pbmc, features = c("Oprd1", "Penk"))

FeaturePlot(pbmc, features = c("Snap25", "Slc17a6", "Slc17a7","Slc6a5"))

## dot plot
features <- c("Snap25", "Oprd1", "Oprm1", "Slc17a6", "Slc17a7", "Gad2")

DotPlot(pbmc, features = features) + RotatedAxis()

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

FeaturePlot(pbmc, features = c("Oprd1", "Penk"))

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

## upset plot of the gene data------
library("UpSetR")

countsData<-read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/gene_data/PG_normalized.csv", 
                     header = T, row.names=1)
gene_interest <- c("Oprd1", "Oprm1", "Penk", "Slc17a6", "Slc17a7")
# "Gad2"

list_input <- vector(mode = "list", length = length(gene_interest))

for (i in seq_along(gene_interest)) {
  gene_name <- countsData %>% 
    .[gene_interest[i], ] 
  
  gene_name1 <- colnames(gene_name)[which(gene_name>0)]
  list_input[[i]] <- gene_name1
  
}

names(list_input) <- gene_interest
upset(fromList(list_input), order.by = "freq")

## analyze gene data from Adam-----
pbmc.data <- read.table("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/gene_data/From\ Adam/cells.count", header = T, sep = "", dec = ".")
pbmc.data <- read.delim("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/gene_data/From\ Adam/cells.count",row.names=1 )

gene_name_data <- read.delim("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/gene_data/From\ Adam/ensembl2entrezID_v90_mm10.txt", )




gene_interest <- c("Oprd1", "Oprm1", "Penk", "Slc17a6", "Slc17a7")
# "Gad2"

list_input <- vector(mode = "list", length = length(gene_interest))

for (i in seq_along(gene_interest)) {
  gene_ID <- which(gene_name_data$Gene.name == gene_interest[i]) %>% 
    gene_name_data$Gene.stable.ID[.]
  
  gene_name <- pbmc.data %>% 
    .[gene_ID, ] 
  
  gene_name1 <- colnames(gene_name)[which(gene_name>0)]
  list_input[[i]] <- gene_name1
  
}

names(list_input) <- gene_interest
upset(fromList(list_input), order.by = "freq")

## for another set of gene
gene_interest <- c( "Slc17a6", "Slc17a7", "Slc32a1","Gad1", "Gad2")
# "Gad2"

list_input <- vector(mode = "list", length = length(gene_interest))

for (i in seq_along(gene_interest)) {
  gene_ID <- which(gene_name_data$Gene.name == gene_interest[i]) %>% 
    gene_name_data$Gene.stable.ID[.]
  
  gene_name <- pbmc.data %>% 
    .[gene_ID, ] 
  
  gene_name1 <- colnames(gene_name)[which(gene_name>0)]
  list_input[[i]] <- gene_name1
  
}

names(list_input) <- gene_interest
upset(fromList(list_input), order.by = "freq")
