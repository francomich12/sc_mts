library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)

#cluster dimplot
DimPlot(data_pt1272, reduction = "umap", label = TRUE)

#dimplot condition|clusters
clusters <- DimPlot(data_pt1272, reduction = 'umap', group.by = 'integrated_snn_res.0.3', label = TRUE)
condition <- DimPlot(data_pt1272, reduction = 'umap', group.by = 'sample')

clusters|condition


#Visualize UMI and gene counts distribution and mitochondrial content
#Visualize UMI and gene counts distribution
VlnPlot(data_pt1272, features = c("nUMI", "nGene"))

VlnPlot(data_pt1272, features = "mitoRatio")

DefaultAssay(data_pt1272) <- "RNA"
DefaultAssay(data_pt1272)

all_markers_pt1272 <-FindAllMarkers(data_pt1272,
                                    logfc.threshold = 0.25,
                                    min.pct = 0.1,
                                    only.pos = TRUE,
                                    test.use = 'wilcox',
                                    slot = 'counts')
# Storing the markers for each cluster in a list
marker_lists_pt1272 <- list()
for(i in 0:7) {   # clusters 0 through 7
  cluster_markers <- all_markers_pt1272 %>% filter(cluster == as.character(i))
  marker_lists_pt1272[[paste0("Cluster ", i)]] <- cluster_markers
  
  # Print the first five markers for each cluster
  cat("Cluster", i, ":", paste(cluster_markers$gene[1:6], collapse = ", "), "\n")
}

#cluster 0
cluster0_fp <- FeaturePlot(data_pt1272, features = c('CAV1', 'TAGLN', 'CAVIN1', 'GNG11', 'MYL9', 'RAMP2'), min.cutoff = 'q10')
cluster0_vp <- VlnPlot(data_pt1272, features = c('CAV1', 'TAGLN', 'CAVIN1', 'GNG11', 'MYL9', 'RAMP2'))
# Cluster 1
FeaturePlot(data_pt1272, features = c('CD24', 'WFDC2', 'AKR1C3', 'MDK', 'CLDN10', 'MARCKSL1'), min.cutoff = 'q10')
VlnPlot(data_pt1272, features = c('CD24', 'WFDC2', 'AKR1C3', 'MDK', 'CLDN10', 'MARCKSL1'))

# Cluster 2
FeaturePlot(data_pt1272, features = c('IL7R', 'CD3D', 'CXCR4', 'CD2', 'TSC22D3', 'BTG1'), min.cutoff = 'q10')
VlnPlot(data_pt1272, features = c('IL7R', 'CD3D', 'CXCR4', 'CD2', 'TSC22D3', 'BTG1'))

# Cluster 3
FeaturePlot(data_pt1272, features = c('CD7', 'NKG7', 'HOPX', 'TNFRSF18', 'TRDC', 'CD247'), min.cutoff = 'q10')
VlnPlot(data_pt1272, features = c('CD7', 'NKG7', 'HOPX', 'TNFRSF18', 'TRDC', 'CD247'))

# Cluster 4
FeaturePlot(data_pt1272, features = c('AIF1', 'C1QC', 'SPI1', 'CD14', 'CD68', 'S100A9'), min.cutoff = 'q10')
VlnPlot(data_pt1272, features = c('AIF1', 'C1QC', 'SPI1', 'CD14', 'CD68', 'S100A9'))

# Cluster 5
FeaturePlot(data_pt1272, features = c('CFAP126', 'TPPP3', 'FAM183A', 'CAPSL', 'C1orf194', 'PIFO'), min.cutoff = 'q10')
VlnPlot(data_pt1272, features = c('CFAP126', 'TPPP3', 'FAM183A', 'CAPSL', 'C1orf194', 'PIFO'))

# Cluster 6
FeaturePlot(data_pt1272, features = c('MALAT1', 'NEAT1', 'PLCG2', 'MUC16', 'IER5L', 'MUC5B'), min.cutoff = 'q10')
VlnPlot(data_pt1272, features = c('MALAT1', 'NEAT1', 'PLCG2', 'MUC16', 'IER5L', 'MUC5B'))

# Cluster 7
FeaturePlot(data_pt1272, features = c('COL6A2', 'LUM', 'DCN', 'MMP2', 'COL6A1', 'CCDC80'), min.cutoff = 'q10')
VlnPlot(data_pt1272, features = c('COL6A2', 'LUM', 'DCN', 'MMP2', 'COL6A1', 'CCDC80'))


cluster_markers <- list(
  c('CAV1', 'TAGLN', 'CAVIN1', 'GNG11', 'MYL9', 'RAMP2'),   # cluster 0
  c('CD24', 'WFDC2', 'AKR1C3', 'MDK', 'CLDN10', 'MARCKSL1'), # Cluster 1
  c('IL7R', 'CD3D', 'CXCR4', 'CD2', 'TSC22D3', 'BTG1'),      # Cluster 2
  c('CD7', 'NKG7', 'HOPX', 'TNFRSF18', 'TRDC', 'CD247'),     # Cluster 3
  c('AIF1', 'C1QC', 'SPI1', 'CD14', 'CD68', 'S100A9'),       # Cluster 4
  c('CFAP126', 'TPPP3', 'FAM183A', 'CAPSL', 'C1orf194', 'PIFO'), # Cluster 5
  c('MALAT1', 'NEAT1', 'PLCG2', 'MUC16', 'IER5L', 'MUC5B'),   # Cluster 6
  c('COL6A2', 'LUM', 'DCN', 'MMP2', 'COL6A1', 'CCDC80')      # Cluster 7
)

#looping through each cluster and generating/ saving plots
for (i in 1:length(cluster_markers)) {
  # Adjust index for cluster number
  cluster_num <- i - 1
  
  #generate FeaturePlot and VlnPlot for each cluster
  fp <- FeaturePlot(data_pt1272, features = cluster_markers[[i]], min.cutoff = 'q10')
  vp <- VlnPlot(data_pt1272, features = cluster_markers[[i]])
  
  #save the plots
  ggsave(paste0("cluster", cluster_num, "_feature_plot.png"), plot = fp, width = 10, height = 8)
  ggsave(paste0("cluster", cluster_num, "_violin_plot.png"), plot = vp, width = 10, height = 8)
}


saveRDS(data_pt1272, file = "/Users/michellefranco/scRNA_patientss/data/pt1272/1272_dge_feature_plots.rds")



