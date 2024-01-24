library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)

#cluster dimplot
DimPlot(data_1315, reduction = "umap", label = TRUE)

#dimplot condition|clusters
clusters <- DimPlot(data_1315, reduction = 'umap', group.by = 'integrated_snn_res.0.3', label = TRUE)
condition <- DimPlot(data_1315, reduction = 'umap', group.by = 'sample')

clusters|condition


#Visualize UMI and gene counts distribution and mitochondrial content
#visualize UMI and gene counts distribution
VlnPlot(data_1315, features = c("nUMI", "nGene"))

VlnPlot(data_1315, features = "mitoRatio")

DefaultAssay(data_1315) <- "RNA"
DefaultAssay(data_1315)

all_markers_pt1315 <-FindAllMarkers(data_1315,
                                    logfc.threshold = 0.25,
                                    min.pct = 0.1,
                                    only.pos = TRUE,
                                    test.use = 'wilcox',
                                    slot = 'counts')
#storing the markers for each cluster in a list
marker_lists_pt1315 <- list()
for(i in 0:7) {   # clusters 0 through 7
  cluster_markers1315 <- all_markers_pt1315 %>% filter(cluster == as.character(i))
  marker_lists_pt1315[[paste0("Cluster1315", i)]] <- cluster_markers1315
  
  #print the first five markers for each cluster
  cat("Cluster1315", i, ":", paste(cluster_markers1315$gene[1:6], collapse = ", "), "\n")
}


#Cluster 0
FeaturePlot(data_1315, features = c('CD24', 'CRABP2', 'CDKN2A', 'KRT18', 'WFDC2', 'EPCAM'), min.cutoff = 'q10')
VlnPlot(data_1315, features = c('CD24', 'CRABP2', 'CDKN2A', 'KRT18', 'WFDC2', 'EPCAM'))
ggsave("cluster0_feature_plot.png", fp0, width = 10, height = 8)
ggsave("cluster0_violin_plot.png", vp0, width = 10, height = 8)

#Cluster 1
FeaturePlot(data_1315, features = c('APOC1', 'RNASE1', 'C1QA', 'C1QB', 'TYROBP', 'AIF1'), min.cutoff = 'q10')
VlnPlot(data_1315, features = c('APOC1', 'RNASE1', 'C1QA', 'C1QB', 'TYROBP', 'AIF1'))
ggsave("cluster1_feature_plot.png", fp1, width = 10, height = 8)
ggsave("cluster1_violin_plot.png", vp1, width = 10, height = 8)

#Cluster 2
FeaturePlot(data_1315, features = c('CD3D', 'CD2', 'LTB', 'TRAC', 'CXCR4', 'CD3E'), min.cutoff = 'q10')
VlnPlot(data_1315, features = c('CD3D', 'CD2', 'LTB', 'TRAC', 'CXCR4', 'CD3E'))
ggsave("cluster2_feature_plot.png", fp2, width = 10, height = 8)
ggsave("cluster2_violin_plot.png", vp2, width = 10, height = 8)

#Cluster 3
FeaturePlot(data_1315, features = c('HOOK2', 'KCNQ1OT1', 'ELF3', 'NUP107', 'IGF2BP2', 'RAB3IP'), min.cutoff = 'q10')
VlnPlot(data_1315, features = c('HOOK2', 'KCNQ1OT1', 'ELF3', 'NUP107', 'IGF2BP2', 'RAB3IP'))
ggsave("cluster3_feature_plot.png", fp3, width = 10, height = 8)
ggsave("cluster3_violin_plot.png", vp3, width = 10, height = 8)

# Cluster 4
FeaturePlot(data_1315, features = c('GNLY', 'CCL5', 'NKG7', 'GZMB', 'CTSW', 'KLRD1'), min.cutoff = 'q10')
VlnPlot(data_1315, features = c('GNLY', 'CCL5', 'NKG7', 'GZMB', 'CTSW', 'KLRD1'))
ggsave("cluster4_feature_plot.png", fp4, width = 10, height = 8)
ggsave("cluster4_violin_plot.png", vp4, width = 10, height = 8)

# Cluster 5
FeaturePlot(data_1315, features = c('ACTA2', 'RGS5', 'CPE', 'PPP1R14A', 'TAGLN', 'MYH11'), min.cutoff = 'q10')
VlnPlot(data_1315, features = c('ACTA2', 'RGS5', 'CPE', 'PPP1R14A', 'TAGLN', 'MYH11'))
ggsave("cluster5_feature_plot.png", fp5, width = 10, height = 8)
ggsave("cluster5_violin_plot.png", vp5, width = 10, height = 8)

# Cluster 6
FeaturePlot(data_1315, features = c('RAMP2', 'ACKR1', 'CAVIN2', 'EMCN', 'MT1X', 'LEPR'), min.cutoff = 'q10')
VlnPlot(data_1315, features = c('RAMP2', 'ACKR1', 'CAVIN2', 'EMCN', 'MT1X', 'LEPR'))
ggsave("cluster6_feature_plot.png", fp6, width = 10, height = 8)
ggsave("cluster6_violin_plot.png", vp6, width = 10, height = 8)

# Cluster 7
FeaturePlot(data_1315, features = c('C1S', 'C1R', 'LUM', 'COL1A2', 'ISLR', 'DCN'), min.cutoff = 'q10')
VlnPlot(data_1315, features = c('C1S', 'C1R', 'LUM', 'COL1A2', 'ISLR', 'DCN'))
ggsave("cluster7_feature_plot.png", fp7, width = 10, height = 8)
ggsave("cluster7_violin_plot.png", vp7, width = 10, height = 8)


cluster_markers1315 <- list(
  c('CD24', 'CRABP2', 'CDKN2A', 'KRT18', 'WFDC2', 'EPCAM'), # Cluster 0
  c('APOC1', 'RNASE1', 'C1QA', 'C1QB', 'TYROBP', 'AIF1'),   # Cluster 1
  c('CD3D', 'CD2', 'LTB', 'TRAC', 'CXCR4', 'CD3E'),         # Cluster 2
  c('HOOK2', 'KCNQ1OT1', 'ELF3', 'NUP107', 'IGF2BP2', 'RAB3IP'), # Cluster 3
  c('GNLY', 'CCL5', 'NKG7', 'GZMB', 'CTSW', 'KLRD1'),       # Cluster 4
  c('ACTA2', 'RGS5', 'CPE', 'PPP1R14A', 'TAGLN', 'MYH11'),   # Cluster 5
  c('RAMP2', 'ACKR1', 'CAVIN2', 'EMCN', 'MT1X', 'LEPR'),     # Cluster 6
  c('C1S', 'C1R', 'LUM', 'COL1A2', 'ISLR', 'DCN')           # Cluster 7
)


for (i in 1:length(cluster_markers1315)) {
  #adjust index for cluster number
  cluster_num1315 <- i - 1
  
  #featurePlot and VlnPlot for each cluster
  fp1315 <- FeaturePlot(data_1315, features = cluster_markers1315[[i]], min.cutoff = 'q10')
  vp1315 <- VlnPlot(data_1315, features = cluster_markers1315[[i]])
  
  #save the plots
  ggsave(paste0("Cluster1315", cluster_num1315, "pt1315_feature_plot.png"), plot = fp, width = 10, height = 8)
  ggsave(paste0("Cluster1315", cluster_num1315, "pt1315_violin_plot.png"), plot = vp, width = 10, height = 8)
}


saveRDS(data_1315, file = "/you_path")



