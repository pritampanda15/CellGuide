#!/usr/bin/env Rscript

# Load required libraries
# Load required libraries
library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(SingleR)
library(celldex)
library(ggplot2)
library(dplyr)

# Command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Please provide the input file path (.h5) and output file path (.rds)")
}

input_file <- args[1]
output_file <- args[2]

# Step 1: Check if the input file is raw 10x Genomics or h5Seurat
if (!grepl(".h5Seurat$", input_file)) {
  message("Detected raw 10x Genomics .h5 file. Reading data and creating a Seurat object...")
  data <- Read10X_h5(input_file)
  seurat_obj <- CreateSeuratObject(counts = data)
} else {
  message("Loading h5Seurat file...")
  seurat_obj <- LoadH5Seurat(input_file)
}

# Step 2: Preprocessing
message("Normalizing data...")
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

message("Identifying variable features...")
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

message("Scaling data...")
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj), verbose = FALSE)

message("Running PCA...")
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), verbose = FALSE)

message("Running UMAP...")
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10, verbose = FALSE)

message("Clustering...")
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Step 3: Save results
message("Saving Seurat object as an RDS file...")
saveRDS(seurat_obj, file = output_file)

# Step 4: Save UMAP plot
message("Generating and saving UMAP plot...")
dir.create("seurat_outputs", showWarnings = FALSE)  # Ensure output directory exists
umap_plot <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("UMAP Plot of Clusters")
ggsave("seurat_outputs/umap_plot.png", umap_plot, width = 8, height = 6)

message("Seurat processing completed.")

# Elbow Plot for PCA
message("Generating and saving PCA elbow plot...")
elbow_plot <- ElbowPlot(seurat_obj, ndims = 50) +
  ggtitle("PCA Elbow Plot")
ggsave("seurat_outputs/pca_elbow_plot.png", elbow_plot, width = 8, height = 6)

# Feature Plot for specific genes
message("Generating and saving feature plots for selected genes...")
genes_of_interest <- c("S100A9", "S100A12", "HBA1", "PPBP", "S100A8", "LYZ")  # Replace with genes relevant to your analysis
for (gene in genes_of_interest) {
  feature_plot <- FeaturePlot(seurat_obj, features = gene) +
    ggtitle(paste("Feature Plot for", gene))
  ggsave(paste0("seurat_outputs/feature_plot_", gene, ".png"), feature_plot, width = 8, height = 6)
}

# Violin Plots for selected genes
message("Generating and saving violin plots for selected genes...")
for (gene in genes_of_interest) {
  violin_plot <- VlnPlot(seurat_obj, features = gene, group.by = "seurat_clusters") +
    ggtitle(paste("Violin Plot for", gene))
  ggsave(paste0("seurat_outputs/violin_plot_", gene, ".png"), violin_plot, width = 8, height = 6)
}

# Heatmap for top markers
# Finding markers and generating heatmap
message("Finding markers and generating heatmap...")
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Group markers by cluster and select the top 10 markers per cluster
top_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Create a heatmap for the top markers
heatmap_plot <- DoHeatmap(seurat_obj, features = top_markers$gene) +
  ggtitle("Heatmap of Top Markers")
ggsave("seurat_outputs/cluster_heatmap.png", heatmap_plot, width = 10, height = 8)

# Dot Plot for selected genes
message("Generating and saving dot plot...")
dot_plot <- DotPlot(seurat_obj, features = genes_of_interest) +
  ggtitle("Dot Plot for Selected Genes")
ggsave("seurat_outputs/dot_plot.png", dot_plot, width = 10, height = 8)

# Run SingleR with updated GetAssayData
message("Running SingleR for cell type annotation...")

# Load reference data
ref_data <- celldex::HumanPrimaryCellAtlasData()

# Run SingleR using the "log-normalized" layer (if you're using Seurat v5 or higher)
singleR_results <- SingleR(
  test = GetAssayData(seurat_obj, layer = "data"), 
  ref = ref_data, 
  labels = ref_data$label.main
)

# Add SingleR labels to Seurat object metadata
seurat_obj[["SingleR.labels"]] <- singleR_results$labels

# Ensure SingleR labels are a factor before assigning to Idents
cell_types <- as.factor(singleR_results$labels)
Idents(seurat_obj) <- cell_types

# Step 3: Save Annotated UMAP Plot
message("Generating and saving UMAP plot with SingleR annotations...")
dir.create("seurat_outputs", showWarnings = FALSE)
umap_plot <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, label.size = 4) +
  ggtitle("UMAP Plot with SingleR Cell Type Annotations")
ggsave("seurat_outputs/umap_singleR_plot.png", umap_plot, width = 8, height = 6)

# Save results
message("Saving Seurat object with SingleR annotations...")
saveRDS(seurat_obj, file = output_file)

message("Annotation with SingleR and plotting completed successfully.")

# Generate Cell Type Statistics
message("Generating cell type statistics...")
cell_type_counts <- as.data.frame(table(singleR_results$labels))
colnames(cell_type_counts) <- c("Cell Type", "Count")
write.csv(cell_type_counts, file = "seurat_outputs/cell_type_summary.csv", row.names = FALSE)

# Generate Differential Gene Analysis Data
message("Identifying differential markers...")
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file = "seurat_outputs/differential_genes.csv", row.names = FALSE)


message("All plots generated and saved successfully.")


