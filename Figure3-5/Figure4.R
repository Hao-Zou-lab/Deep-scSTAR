#--------------------------------------------------------------
# filename : Figure4.R
# Date : 2024-8-29
# contributor : Lianchong Gao
# Python version: 4.2.2
# Platform: Windows PC
#--------------------------------------------------------------

library(CellChat)
library(Seurat)
library(tidyverse)
library(viridis)
library(RColorBrewer)

setwd('~/Fig4S')
Sobj <- readRDS('Fig4_Sobj.rds')
DefaultAssay(Sobj) <- "SCT"

###############################################################################
#'                             Manuscript: Figure 4d                   '#
###############################################################################
cluster_colors <- c(
  "C_0" = "#bbdf9e",  
  "C_1" = "#3f7fb5",   
  "C_2" = "#5baeab",   
  "C_3" = "#f1945a",   
  "C_4" = "#f1dc92",   
  "C_5" = "#9d144c",   
  "C_6" = "#60539e",   
  "C_7" = "#f5bc70"     
)

plot <- SpatialDimPlot(Sobj,
                       group.by = "Dsc_Cluster",
                       pt.size.factor = 6.5,
                       image.alpha=0.1,
                       cols = cluster_colors, 
                       #stroke =NA, 
                       images = "imga15") +
  theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 16)))

print(plot)
ggsave("Dsc_SpatialDim.png", plot = plot, dpi = 600, width =5, height = 5, units = "in")

cluster_colors <- c(
  "C_0" = "#bbdf9e",  
  "C_1" = "#3f7fb5",   
  "C_2" = "#5baeab",   
  "C_3" = "#f1945a",   
  "C_4" = "#f1dc92",   
  "C_5" = "#9d144c",   
  "C_6" = "#60539e",   
  "C_7" = "#f5bc70"     
)

plot <- SpatialDimPlot(Sobj,
                       group.by = "Raw_Cluster",
                       pt.size.factor = 6.5,
                       image.alpha=0.1,
                       cols = cluster_colors, 
                       #stroke =NA, 
                       images = "imga15") +
  theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 16)))
print(plot)
ggsave("Raw_SpatialDim.png", plot = plot, dpi = 600, width =5, height = 5, units = "in")

###############################################################################
#'                             Manuscript: Figure 4e                   '#
###############################################################################
gene_sets_list <- c('COL1A1','COL3A1','ACTA2','RGS5','MCAM','PDGFRB','GZMK','CD3D','TCF7','IL7R','FOSB','AEBP1')
cluster_levels <- c("C_0", "C_1", "C_5", "C_3", "C_2", "C_4", "C_6", "C_7")
Sobj$Dsc_Cluster <- factor(Sobj$Dsc_Cluster, levels = cluster_levels)
dotplot <- DotPlot(Sobj, features = gene_sets_list, group.by = "Dsc_Cluster") +
  scale_colour_gradientn(colors = c("#4897c5", "#f3f3f3", "#c4383e")) +  
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size = 6),
        axis.text.y = element_text(size = 6),
        legend.position="1") +
  ggtitle("")
print(dotplot)
ggsave("Dsc_DotMarker.png", plot = dotplot, dpi = 600, width =5, height = 5, units = "in")

###############################################################################
#'                             Manuscript: Figure 4f                   '#
###############################################################################
gene_sets_list <- list('M1_TAM' = c('IL1B','TNF','CXCL9','CXCL10','CXCL11','FCGR1A','IRF1','HLA-DPB1','CD86','MARCO','IL2RA'),
                       'M2_TAM' = c('CXCR4','IL27RA','CSF1R','CCL7','CCL17','CCL23','CTSD','GAS7','HMOX1',
                                'PPARG','LIPA','CD209','CLEC7A','F13A1','MAF','MS4A4A'),
                       'Progenitor Exhausted CD8' = c("IL7R", "GPR183", "LMNA", "NR4A3", "TCF7", "MGAT4A", 
                                                      "CD55", "AIM1", "PER1", "FOSL2", "EGR1", "TSPYL2", "YPEL5", 
                                                      "CSRNP1", "REL", "SKIL", "PIK3R1", "FOXP1", "RGCC", "PFKFB3", 
                                                      "MYADM", "ZFP36L2", "USP36", "TC2N", "FAM177A1", "BTG2", "TSC22D2", 
                                                      "FAM65B", "STAT4", "RGPD5", "NEU1", "IFRD1", "PDE4B", "NR4A1"),
                       'Terminally Exhausted CD8' = c(
                         "NKG7", "RAC2", "CLIC1", "GZMA", "PRF1", "APOBEC3C", "RHOA", "CCL4", "COTL1", "PSME2",
                         "HLA-DPA1", "HMGN2", "LSP1", "PSMB9", "LCK", "SRP14", "ARPC3", "ARPC1B", "TPI1", "APOBEC3G",
                         "HLA-DPB1", "LDHB", "ATP5G2", "MYL12B", "PSMB8", "PSMA7", "HLA-DRB1", "SUB1", "ARPC4", "CTSW",
                         "SUMO2", "TAP1", "GZMB", "RARRES3", "CAP1", "UCP2", "PPIB", "RAN", "CHCHD2", "PARK7",
                         "HCST", "GABARAP", "HLA-DRA", "SOD1", "CAPZB", "S100A4", "RNASEK", "PPP1CA", "PKM", "IFI16",
                         "ACTR3", "ITM2A", "SLC25A5", "PGAM1", "ANXA6", "CD27", "ATP5B", "LYST", "PSMB10", "MIF",
                         "LY6E", "ANKRD10", "CTSD", "UBE2L6", "EDF1", "NONO", "TIGIT", "FKBP1A", "IL2RB", "HMGN1",
                         "ATP5L", "GZMH", "STAT1", "GPI", "LCP2", "GBP2", "ARL6IP5", "CCL4L1", "PRDM1", "OST4",
                         "PDCD1", "HINT1", "HNRNPF", "GBP5", "COX7C", "ARPC5", "GIMAP4", "XRCC6", "C17orf62", "PRR13",
                         "HLA-DRB5", "WDR1", "ARL6IP1", "ISG15", "ATP5A1", "EWSR1", "COPE", "HAVCR2", "EIF3H", "ANXA5",
                         "C11orf58", "IFI6", "SIRPG", "CALM3", "SHISA5", "DENND2D", "MAP4K1", "BUB3", "IKZF3", "SNRPB",
                         "EID1", "PSMB1", "PTPN6", "NDUFA13", "SSR4", "COX8A", "PTPN7", "MAT2B", "PSTPIP1", "GSTP1",
                         "PSMB3", "IRF9", "TRAF3IP3", "GIMAP7", "PSMA2", "SASH3", "CD164", "ETNK1", "S100A11", "KLRD1",
                         "MOB1A", "SH2D1A", "UBE2V1", "SH3KBP1", "ATP5G3", "PSMA5", "MT2A", "LAT", "IFNG", "RAB27A",
                         "COX5A", "DDOST", "PSMB4", "SRP9", "BRK1", "TNFRSF1B", "EIF4H", "GMFG", "ANXA2", "TCEB2",
                         "RBPJ", "COX6A1", "UBXN1", "PSMD8", "CD63", "ATP6V0E1", "NDUFB8", "CTSC", "SNRPD2", "ATP5C1",
                         "PRELID1", "COX7A2", "PSMA6", "ECH1", "U2AF1", "HMGB2", "FAM49B", "CD38", "TSPO", "IDH2",
                         "CASP4", "CCL3", "TRMT112", "SURF4", "PSMA1", "YWHAE", "LASP1", "PYHIN1", "ANAPC16", "TUBB",
                         "CSNK2B", "PRKAR1A", "SLAMF7", "GNG5", "PRDX1", "RQCD1", "CCNDBP1", "INPP4B", "CDK2AP2", "CBX3",
                         "RPN1", "SPCS1", "PSMA3", "SIT1", "XRCC5", "EIF3C", "CXCR6", "COX6C", "M6PR", "ANP32E"
                       ))
Sobj <- SCTransform(Sobj, assay = "Spatial", verbose = T)
expr_matrix <- GetAssayData(Sobj, assay = "SCT", slot = "data")
expr_matrix <- as.matrix(expr_matrix)
gsva_scores <- gsva(expr_matrix, gene_sets_list, kcdf="Gaussian",method = "gsva", verbose = FALSE)
gsva_scores_df <- as.data.frame(t(gsva_scores))
rownames(gsva_scores_df) <- colnames(expr_matrix)
Sobj <- AddMetaData(Sobj, metadata = gsva_scores_df)


Obj1 <- subset(Sobj, subset = Dsc_Cluster %in% c('C_1', 'C_5', 'C_3'))
cluster_levels <- c("C_1", "C_5", "C_3")  
Obj1$Dsc_Cluster <- factor(Obj1$Dsc_Cluster, levels = cluster_levels)
data_to_plot <- FetchData(Obj1, vars = c("Dsc_Cluster", "Progenitor Exhausted CD8", "Terminally Exhausted CD8", "M1_TAM", "M2_TAM"))
data_long <- reshape2::melt(data_to_plot, id.vars = "Dsc_Cluster")
plot <- ggplot(data_long, aes(x = Dsc_Cluster, y = value, fill = Dsc_Cluster)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  facet_wrap(~ variable, scales = "free_y") +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  scale_fill_manual(values = c("#5d9ace", "#b63469", "#ffaa75")) +
  labs(title = "", x = "Dsc_Cluster", y = "Expression Score") +
  theme_classic() +
  theme(text = element_text(color = "black", face = "bold"),  
        axis.title = element_text(size = 12,color = "black", face = "bold"), 
        axis.text = element_text(size = 10,color = "black", face = "bold"),  
        strip.text = element_text(size = 10,color = "black", face = "bold"),
        #legend.position="top"
  )  

print(plot)
ggsave("Figure4f.png", plot = plot, dpi = 600, width =8, height = 4.5, units = "in")

###############################################################################
#'                             Manuscript: Figure 4g                   '#
###############################################################################
R2 <- subset(Sobj, subset = Dsc_Cluster %in% c('C_2', 'C_4', 'C_6'))
gene_sets <- list(
  MSGs = c('MCAM','PDGFA','MEF2C','CSPG4'),
  Tcell_TAM = c('GZMK','DOCK8','CD14','S100A12','CD8A'),
  Plasma = c('JCHAIN','PIGR','IGKC')
)
dotplot <- DotPlot(R2, features = gene_sets, group.by = "Dsc_Cluster") +
  scale_colour_gradientn(colors = c("#4897c5", "#f3f3f3", "#c4383e")) +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = "right") +
  ggtitle("Gene Expression DotPlot by Cluster")
print(dotplot)
ggsave("Dsc_DotMarker2.png", plot = dotplot, dpi = 600, width =5, height = 5, units = "in")

###############################################################################
#'                             Manuscript: Figure 4h&i                   '#
###############################################################################
Idents(Sobj) <- "Dsc_Cluster"
head(Sobj)

color.use <- scPalette(nlevels(Sobj))
names(color.use) <- levels(Sobj)

data.input = Seurat::GetAssayData(Sobj, slot = "data", assay = "SCT") 
meta = data.frame(labels = Idents(Sobj), 
                  row.names = names(Idents(Sobj))) 
meta$labels <- factor(meta$labels,
                      levels = c("C_0", "C_1", "C_5", "C_3", "C_2", "C_4", "C_6", "C_7"))

unique(meta$labels)
spatial.locs = Seurat::GetTissueCoordinates(Sobj, scale = NULL, 
                                            cols = c("imagerow", "imagecol")) 
existing_scale_factors <- Sobj@images$imga15@scale.factors

scale.factors <- list(
  spot.diameter = existing_scale_factors$spot, 
  spot = existing_scale_factors$spot, 
  fiducial = existing_scale_factors$fiducial,
  hires = existing_scale_factors$hires,
  lowres = existing_scale_factors$lowres
)
#print(scale.factors)

cellchat <- createCellChat(object = data.input, 
                           meta = meta, 
                           group.by = "labels", 
                           datatype = "spatial", ###
                           coordinates = spatial.locs, 
                           scale.factors = scale.factors)




CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) 
future::plan("multisession", workers = 8) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)


cellchat <- computeCommunProb(cellchat, 
                              type = "truncatedMean", trim = 0.1, 
                              distance.use = TRUE, 
                              scale.distance = 0.01)
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)



groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions",color.use =cluster_colors)
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use =cluster_colors)

cellchat@netP$pathways
pathways.show  <- 'CD99'
vertex.receiver = c(6,4,3,5)  
netVisual_aggregate(cellchat, signaling = pathways.show,                      
                    vertex.receiver = vertex.receiver,layout = "hierarchy")

pathways.show  <- 'FN1'
vertex.receiver = c(6,4,3,5) 
netVisual_aggregate(cellchat, signaling = pathways.show,                      
                    vertex.receiver = vertex.receiver,layout = "hierarchy")


