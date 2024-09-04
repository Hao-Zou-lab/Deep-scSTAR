#--------------------------------------------------------------
# filename : Figure5.R
# Date : 2024-8-30
# contributor : Lianchong Gao
# Python version: 4.2.2
# Platform: Windows PC
#--------------------------------------------------------------
setwd("~/Fig5/")
library(pls)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(MCPcounter)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(devtools)
library(ROCit)
library(ggsignif)
library(spacexr)
library(ggtext)
library(tidyr)
library(SCpubr)
library(circlize)
library(utils)
library(rjson)
library(limma)
library(estimate)
library(GSVA)
library(parallel)
library(readxl)
library(scales) 
library(ggpubr)
library(ggExtra)

numCores <- parallel::detectCores()
print(numCores)
###############################################################################
#'                             Manuscript: Figure 5a&b                   '#
###############################################################################

HCC_Dsc <- readRDS('HCC_Dsc.rds')
HCC_Ori <- readRDS('HCC_Ori.rds')

colors <- c("Neu_C0" = "#5ec4af","Neu_C1" = "#fdc785","Neu_C2" = "#bb405d","Neu_C3" = "#2d3293","Neu_C4" = "#f4cf66","Neu_C5" = "#bce2ee","Neu_C6" = "#367ba2")
p <- SCpubr::do_DimPlot(HCC_Dsc,colors.use = colors,group.by = "Dsc_cluster",legend.position = "right")
ggsave("UMAP_Dsc.png", plot = p, width = 6, height = 5, dpi = 600)

colors <- c("Raw_C0" = "#5ec4af","Raw_C1" = "#fdc785","Raw_C2" = "#bb405d","Raw_C3" = "#2d3293","Raw_C4" = "#f4cf66","Raw_C5" = "#bce2ee","Raw_C6" = "#367ba2")
p <- SCpubr::do_DimPlot(HCC_Ori,colors.use = colors,group.by = "Raw_cluster",legend.position = "right")
ggsave("UMAP_Raw.png", plot = p, width = 6, height = 5, dpi = 600)

rm(HCC_Dsc)

genes <- list("VEGFA+ Neu" = c("VEGFA","SLC3A2", "P4HA1","MIF","RALGDS"),
              "CXCR2+ Neu" = c("CXCR2","MOSPD2","MRVI1","UBR2","CXCR1"),
              "S100A12+ Neu" = c("S100A12","S100A4","CDA","S100A6",'S100A8'))

custom_order_dsc <- c("Neu_C3","Neu_C2","Neu_C4", "Neu_C1", "Neu_C5", "Neu_C6",  "Neu_C0") 
custom_order_raw <- c("Raw_C3", "Raw_C1", "Raw_C2", "Raw_C6", "Raw_C4", "Raw_C5", "Raw_C0") 

HCC_Ori$Dsc_cluster <- factor(HCC_Ori$Dsc_cluster, levels = custom_order_dsc)
HCC_Ori$Raw_cluster <- factor(HCC_Ori$Raw_cluster, levels = custom_order_raw)

Idents(HCC_Ori) <- HCC_Ori$Dsc_cluster
p1 <-SCpubr::do_DotPlot(HCC_Ori, 
                        features = genes, 
                        #cluster = TRUE,
                        dot.scale = 10,
                        cluster.idents = FALSE,
                        colors.use = c("#6fa7ce",'#f4f8f9', "#c4383e"))
p1
ggsave("Dot_Dsc.png", plot = p1, width = 6, height = 5, dpi = 600)

Idents(HCC_Ori) <- HCC_Ori$Raw_cluster
p2 <-SCpubr::do_DotPlot(HCC_Ori, 
                        features = genes, 
                        #cluster = TRUE,
                        cluster.idents = FALSE,
                        dot.scale = 10,
                        colors.use = c("#6fa7ce",'#f4f8f9', "#c4383e"))
p2
ggsave("Dot_Raw.png", plot = p2, width = 6, height = 5, dpi = 600)


###############################################################################
#'                             Manuscript: Figure 5c                   '#
###############################################################################
table0 <- table(HCC_Ori$Dsc_cluster,HCC_Ori$Raw_cluster)
M <- colSums(table0)
N <- rowSums(table0)
K <- length(HCC_Ori$Raw_cluster)
P <- matrix(1, length(N), length(M))
for (i in 1:length(N)){
  for ( j in 1:length(M)){
    A <- matrix(0,2,2)
    A[1,1] <- table0[i,j]
    A[1,2] <- N[i]-table0[i,j]
    A[2,1] <- M[j]-table0[i,j]
    A[2,2] <- N[i]+M[j]-table0[i,j]
    Ar <- rowSums(A)
    Ac <- colSums(A)
    P[i, j] <- 1-phyper(table0[i,j], M[j], K-M[j], N[i], lower.tail = TRUE, log.p = FALSE)
  }
}
colnames(P) <- colnames(table0)
rownames(P) <- rownames(table0)
my.data = table0

colors <- c("0" = "#5ec4af",
            "1" = "#fdc785",
            "2" = "#bb405d",
            "3" = "#2d3293",
            "4" = "#f4cf66",
            "5" = "#bce2ee",
            "6" = "#367ba2")

col1 <- setNames(colors, paste("Neu_C", 0:6, sep=""))
col2 <- setNames(colors, paste("Raw_C", 0:6, sep=""))

coul <- col1
grid.col=NULL
grid.col[colnames(my.data)]=col2
grid.col[rownames(my.data)]=coul
circos.clear()
circos.par(gap.degree=c(rep(2,nrow(my.data)-1),10,rep(2,ncol(my.data)-1),10),start.degree=180)
my.P = P

linkc1 = col2
linkc = matrix("white",dim(my.P)[1], dim(my.P)[2])
for (j in 1:dim(my.P)[2]){
  for (i in 1:dim(my.P)[1]){
    if (my.P[i,j]<0.05){
      linkc[i,j] = linkc1[j]
    }
  }
}

##########plot
png("Cirplot.png", width=550.5, height=550.5)
par(cex = 1.5, mar = c(0, 0, 0, 0))
# 绘制弦图
chordDiagram(my.data, directional = TRUE, grid.col = grid.col, col = linkc,
             transparency = 0.8, diffHeight = 0.06,
             annotationTrack = c("name", "grid"),
             annotationTrackHeight = c(0.05, 0.05))
dev.off()

###############################################################################
#'                             Manuscript: Figure 5e                   '#
###############################################################################
DscMarker <- read_excel("~/HCC/TCGA-LIHC/Neutrophil_Markersold.xlsx") # Dsc Marker genes
RawMarker <- read.csv("~/scRNAseqData of HCC/Neutrophil_Markers.csv") # Raw Marker genes

DscMarker <- DscMarker %>%
  filter(`Adjusted.P.value` < 0.01) %>%
  arrange(desc(`Log2.Fold.Change`))
top_genes_by_cluster <- DscMarker %>%
  group_by(Cluster) %>%
  slice_head(n = 20) %>%
  summarise(Genes = list(Gene), .groups = 'drop')
top_genes_by_cluster <- top_genes_by_cluster %>%
  mutate(Cluster = paste(Cluster))
DscMarker_list <- setNames(top_genes_by_cluster$Genes, top_genes_by_cluster$Cluster)

RawMarker <- RawMarker %>%
  filter(`Adjusted.P.value` < 0.01) %>%
  arrange(desc(`Log2.Fold.Change`))
top_genes_by_cluster <- RawMarker %>%
  group_by(Cluster) %>%
  slice_head(n = 20) %>%
  summarise(Genes = list(Gene), .groups = 'drop')
top_genes_by_cluster <- top_genes_by_cluster %>%
  mutate(Cluster = paste(Cluster))
RawMarker_list <- setNames(top_genes_by_cluster$Genes, top_genes_by_cluster$Cluster)


# read RNAseq data
file_path <- "~/HCC/GSE140901_processed_data.txt"
data <- read.table(file_path, header = TRUE, sep = "\t", comment.char = "", check.names = FALSE, row.names = 1, skip = 4)

data <- data[, colSums(is.na(data)) != nrow(data)]

sample_status <- rep("NR", ncol(data))
known_R_samples <- c("S30", "S22", "S19", "S18")
sample_status[colnames(data) %in% known_R_samples] <- "R"

meta_data <- data.frame(
  sample = colnames(data),
  Group = sample_status,
  row.names = colnames(data)
)
data_matrix <- as.matrix(data)

# gsva
gsva_results <- gsva(data_matrix, gset.idx.list=DscMarker_list, 
                     kcdf="Gaussian", #"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
                     verbose=TRUE, 
                     parallel.sz =  parallel::detectCores())

df_gsva <- as.data.frame(t(gsva_results))
colnames(df_gsva) <- rownames(gsva_results)
df_gsva$sample <- rownames(df_gsva)

df_plot <- merge(df_gsva, meta_data, by = "sample")
df_plot_selected <- df_plot[, c("sample", "Group", "Neu_C_1")]
colnames(df_plot_selected)[3] <- "value"

shapiro.test(df_plot_selected$value[df_plot_selected$Group == "R"])
shapiro.test(df_plot_selected$value[df_plot_selected$Group == "NR"])
var.test(value ~ Group, data = df_plot_selected)

p_value <- t.test(value ~ Group, data = df_plot_selected, var.equal = FALSE)$p.value
p_value

df_plot_long <- pivot_longer(df_plot_selected, cols = "value", names_to = "gene_set", values_to = "value")

p <- ggplot(df_plot_long, aes(x = Group, y = value, color = Group, fill = Group)) +
  geom_boxplot(fill = NA, size = 1.2) +  
  geom_jitter(position = position_jitter(width = 0.2), aes(color = Group), size = 1) + 
  labs(title = "",
       x = "response",
       y = "GSVA Score of Neu_C1") +
  theme_minimal() +
  scale_fill_manual(values = c("NR" = "#cd5150", "R" = "#81b1d3")) +  
  scale_color_manual(values = c("NR" = "#cd5150", "R" = "#81b1d3")) +  
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),  
    axis.text = element_text(color = "black"),  
    axis.title = element_text(color = "black"),  
    axis.text.x = element_text(angle = 45, hjust = 1)  
  ) 
print(p)

#Figure 5e   
ggsave(filename = "gsva_NeuC1_NRvsR.png", plot = p, dpi = 600, width = 3, height = 3)

##########################################################################################################
#FigS11

# gsva
gsva_results <- gsva(data_matrix, gset.idx.list=RawMarker_list, 
                     kcdf="Gaussian", #"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
                     verbose=TRUE, 
                     parallel.sz =  parallel::detectCores())

df_gsva <- as.data.frame(t(gsva_results))
colnames(df_gsva) <- rownames(gsva_results)
df_gsva$sample <- rownames(df_gsva)

df_plot <- merge(df_gsva, meta_data, by = "sample")
df_plot_selected <- df_plot[, c("sample", "Group", "Raw_C_2")]
colnames(df_plot_selected)[3] <- "value"

p_value <- t.test(value ~ Group, data = df_plot_selected, var.equal = FALSE)$p.value
p_value

df_plot_long <- pivot_longer(df_plot_selected, cols = "value", names_to = "gene_set", values_to = "value")

p <- ggplot(df_plot_long, aes(x = Group, y = value, color = Group, fill = Group)) +
  geom_boxplot(fill = NA, size = 1.2) +  # 增粗边框线
  geom_jitter(position = position_jitter(width = 0.2), aes(color = Group), size = 1) + 
  labs(title = "",
       x = "response",
       y = "GSVA Score of Raw_C2") +
  theme_minimal() +
  scale_fill_manual(values = c("NR" = "#cd5150", "R" = "#81b1d3")) +  
  scale_color_manual(values = c("NR" = "#cd5150", "R" = "#81b1d3")) +  
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),  
    axis.text = element_text(color = "black"),  
    axis.title = element_text(color = "black"),  
    axis.text.x = element_text(angle = 45, hjust = 1)  
  ) 
print(p)
ggsave(filename = "gsva_RawC2_NRvsR.png", plot = p, dpi = 600, width = 3, height = 3)

df_plot_selected <- df_plot[, c("sample", "Group", "Raw_C_6")]
colnames(df_plot_selected)[3] <- "value"
p_value <- t.test(value ~ Group, data = df_plot_selected, var.equal = FALSE)$p.value
p_value

df_plot_long <- pivot_longer(df_plot_selected, cols = "value", names_to = "gene_set", values_to = "value")
p <- ggplot(df_plot_long, aes(x = Group, y = value, color = Group, fill = Group)) +
  geom_boxplot(fill = NA, size = 1.2) +  
  geom_jitter(position = position_jitter(width = 0.2), aes(color = Group), size = 1) + 
  labs(title = "",
       x = "response",
       y = "GSVA Score of Raw_C6") +
  theme_minimal() +
  scale_fill_manual(values = c("NR" = "#cd5150", "R" = "#81b1d3")) +  
  scale_color_manual(values = c("NR" = "#cd5150", "R" = "#81b1d3")) +  
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),  
    axis.text = element_text(color = "black"),  
    axis.title = element_text(color = "black"),  
    axis.text.x = element_text(angle = 45, hjust = 1)  
  ) 
print(p)
ggsave(filename = "gsva_RawC6_NRvsR.png", plot = p, dpi = 600, width = 3, height = 3)





shapiro.test(df_plot_selected$value[df_plot_selected$Group == "R"])
shapiro.test(df_plot_selected$value[df_plot_selected$Group == "NR"])
var.test(value ~ Group, data = df_plot_selected)


###############################################################################
#'                             Manuscript: Figure 5f                   '#
###############################################################################

metafile="D:\\ccrcc\\HCC\\TCGA-LIHC\\metadata.cart.2024-07-30T.json" 
gdcfliename="D:\\ccrcc\\HCC\\TCGA-LIHC\\gdc_download_20240730_093629.670989" 
path1="D:\\ccrcc\\HCC\\TCGA-LIHC\\gdc_download_20240730_093629.670989\\" 
outfilename="TCGA-LIHC_FPKM_Tumor.txt" #

json = jsonlite::fromJSON(metafile)
id = json$associated_entities[[1]][,1]
sample_id = sapply(json$associated_entities,function(x){x[,1]})
file_sample = data.frame(sample_id,file_name=json$file_name)  
count_file <- list.files(gdcfliename,pattern = '*gene_counts.tsv',recursive = TRUE)
count_file_name <- strsplit(count_file,split='/')
count_file_name <- sapply(count_file_name,function(x){x[2]})
matrix = data.frame(matrix(nrow=60660,ncol=0))
for (i in 1:length(count_file_name)){
  path = paste0(path1,count_file[i])
  data<- read.delim(path,fill = TRUE,header = FALSE,row.names = 1)
  colnames(data)<-data[2,]
  data <-data[-c(1:6),]
  data <- data[7] 
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix <- cbind(matrix,data)
}
sample1 = paste0(path1,count_file[1])
names=read.delim(sample1,fill = TRUE,header = FALSE,row.names = 1)
colnames(names)<-names[2,]
names <-names[-c(1:6),]
names = names[,1:2]
same=intersect(rownames(matrix),rownames(names))
matrix=matrix[same,]
names=names[same,]
matrix$symbol=names[,1]
matrix=matrix[,c(ncol(matrix),1:(ncol(matrix)-1))]
write.table(matrix,file=outfilename,row.names = F,quote = F,sep = "\t")

matrix_Normal <- read.table("D:\\ccrcc\\HCC\\TCGA-LIHC\\TCGA-LIHC_FPKM_Normal.txt", header=TRUE, sep="\t")
non_duplicate_rows <- !duplicated(matrix_Normal[, 1])
matrix_Normal <- matrix_Normal[non_duplicate_rows, ]
rownames(matrix_Normal) <- matrix_Normal[, 1]
matrix_Normal <- matrix_Normal[, -1]
head(matrix_Normal)

matrix_Tumor <- read.table("D:\\ccrcc\\HCC\\TCGA-LIHC\\TCGA-LIHC_FPKM_Tumor.txt", header=TRUE, sep="\t")
non_duplicate_rows <- !duplicated(matrix_Tumor[, 1])
matrix_Tumor <- matrix_Tumor[non_duplicate_rows, ]
rownames(matrix_Tumor) <- matrix_Tumor[, 1]
matrix_Tumor <- matrix_Tumor[, -1]
head(matrix_Tumor)

sample_names_Normal <- colnames(matrix_Normal)
sample_names_Tumor <- colnames(matrix_Tumor)
tissue_type_Normal <- rep("N", length(sample_names_Normal))
tissue_type_Tumor <- rep("T", length(sample_names_Tumor))
meta.data_Normal <- data.frame(sample = sample_names_Normal, Tissue = tissue_type_Normal)
meta.data_Tumor <- data.frame(sample = sample_names_Tumor, Tissue = tissue_type_Tumor)
meta.data <- rbind(meta.data_Normal, meta.data_Tumor)
head(meta.data)



input.f <- "D:\\ccrcc\\HCC\\TCGA-LIHC\\TCGA-LIHC_FPKM_Tumor_Cleaned.txt"
filterCommonGenes(input.f=input.f, output.f="estimate_input.gct", id="GeneSymbol")#

estimateScore("estimate_input.gct", "estimate_score.gct", platform="affymetrix")#
plotPurity(scores="estimate_score.gct", platform="affymetrix", output.dir = "figures/")
gct_data <- read.table("estimate_score.gct", header=TRUE, sep="\t", comment.char="#")
head(gct_data)
gct_data <- read.table("estimate_score.gct", header=TRUE, sep="\t", check.names=FALSE, comment.char="#")

sample_names <- gct_data[1, -c(1, 2)]  
sample_names <- unlist(sample_names)  

names(sample_names) <- seq_along(sample_names)
StromalScore <- as.numeric(gct_data[2, -c(1, 2)]) 

mean_StromalScore <- mean(StromalScore)
group_labels <- ifelse(StromalScore > mean_StromalScore, "high", "low")
meta.data <- data.frame(sample = sample_names, Group = group_labels)
print(meta.data)


common_rows <- intersect(rownames(matrix_Normal), rownames(matrix_Tumor))
matrix_Normal_common <- matrix_Normal[common_rows, ]
matrix_Tumor_common <- matrix_Tumor[common_rows, ]
combined_matrix <- cbind(matrix_Normal_common, matrix_Tumor_common)
matrix_Tumor <- as.matrix(matrix_Tumor)
head(combined_matrix)

combined_matrix <- as.matrix(combined_matrix)

gsva_results <- gsva(matrix_Tumor, gset.idx.list=DscMarker_list, 
                     kcdf="Gaussian", #"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
                     verbose=T, 
                     parallel.sz = parallel::detectCores())

df_gsva <- as.data.frame(t(gsva_results))
colnames(df_gsva) <- rownames(gsva_results)
df_gsva$sample <- rownames(df_gsva)

df_plot <- merge(df_gsva, meta.data, by = "sample")
df_plot_selected <- df_plot[, c("sample", "Group", "Neu_C_1")]

shapiro.test(df_plot_selected$Neu_C_1[df_plot_selected$Group == "low"])
shapiro.test(df_plot_selected$Neu_C_1[df_plot_selected$Group == "high"])
var.test(Neu_C_1 ~ Group, data = df_plot_selected)

library(tidyr)
df_plot_long <- pivot_longer(df_plot_selected, cols = "Neu_C_1")

wilcox_test_result <- wilcox.test(Neu_C_1 ~ Group, data = df_plot_selected)
print(wilcox_test_result)

p <- ggplot(df_plot_long, aes(x = Group, y = value, color = Group, fill = Group)) +
  geom_boxplot(fill = NA, size = 1.2) +  
  geom_jitter(position = position_jitter(width = 0.2), aes(color = Group), size = 1) +  
  labs(title = "",
       x = "StromalScore",
       y = "GSVA Score of Neu_C1") +
  theme_minimal() +
  scale_fill_manual(values = c("high" = "#cd5150", "low" = "#81b1d3")) +  
  scale_color_manual(values = c("high" = "#cd5150", "low" = "#81b1d3")) +  
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),  
    axis.text = element_text(color = "black"),  
    axis.title = element_text(color = "black"),  
    axis.text.x = element_text(angle = 45, hjust = 1) 
  ) 

print(p)
#Figure 5f  
ggsave(filename = "gsva_StromalScore_NeuC1.png", plot = p, dpi = 600, width = 3, height = 3)

#################################################################################################
#FigS11
gsva_results <- gsva(matrix_Tumor, gset.idx.list=RawMarker_list, 
                     kcdf="Gaussian", #"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
                     verbose=T, 
                     parallel.sz = parallel::detectCores())

df_gsva <- as.data.frame(t(gsva_results))
colnames(df_gsva) <- rownames(gsva_results)
df_gsva$sample <- rownames(df_gsva)

df_plot <- merge(df_gsva, meta.data, by = "sample")
df_plot_selected <- df_plot[, c("sample", "Group", "Raw_C_2")]

shapiro.test(df_plot_selected$Raw_C_2[df_plot_selected$Group == "low"])
shapiro.test(df_plot_selected$Raw_C_2[df_plot_selected$Group == "high"])
var.test(Raw_C_2 ~ Group, data = df_plot_selected)

library(tidyr)
df_plot_long <- pivot_longer(df_plot_selected, cols = "Raw_C_2")

wilcox_test_result <- wilcox.test(Raw_C_2 ~ Group, data = df_plot_selected)
print(wilcox_test_result)

p <- ggplot(df_plot_long, aes(x = Group, y = value, color = Group, fill = Group)) +
  geom_boxplot(fill = NA, size = 1.2) +  
  geom_jitter(position = position_jitter(width = 0.2), aes(color = Group), size = 1) +  
  labs(title = "",
       x = "StromalScore",
       y = "GSVA Score of Raw_C_2") +
  theme_minimal() +
  scale_fill_manual(values = c("high" = "#cd5150", "low" = "#81b1d3")) +  
  scale_color_manual(values = c("high" = "#cd5150", "low" = "#81b1d3")) +  
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),  
    axis.text = element_text(color = "black"),  
    axis.title = element_text(color = "black"),  
    axis.text.x = element_text(angle = 45, hjust = 1) 
  ) 

print(p)
ggsave(filename = "gsva_StromalScore_RawC2.png", plot = p, dpi = 600, width = 3, height = 3)



df_plot_selected <- df_plot[, c("sample", "Group", "Raw_C_6")]

shapiro.test(df_plot_selected$Raw_C_6[df_plot_selected$Group == "low"])
shapiro.test(df_plot_selected$Raw_C_6[df_plot_selected$Group == "high"])
var.test(Raw_C_6 ~ Group, data = df_plot_selected)

library(tidyr)
df_plot_long <- pivot_longer(df_plot_selected, cols = "Raw_C_6")

wilcox_test_result <- wilcox.test(Raw_C_6 ~ Group, data = df_plot_selected)
print(wilcox_test_result)

p <- ggplot(df_plot_long, aes(x = Group, y = value, color = Group, fill = Group)) +
  geom_boxplot(fill = NA, size = 1.2) +  
  geom_jitter(position = position_jitter(width = 0.2), aes(color = Group), size = 1) + 
  labs(title = "",
       x = "StromalScore",
       y = "GSVA Score of Raw_C_6") +
  theme_minimal() +
  scale_fill_manual(values = c("high" = "#cd5150", "low" = "#81b1d3")) +  
  scale_color_manual(values = c("high" = "#cd5150", "low" = "#81b1d3")) +  
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),  
    axis.text = element_text(color = "black"),  
    axis.title = element_text(color = "black"),  
    axis.text.x = element_text(angle = 45, hjust = 1) 
  ) 

print(p)
ggsave(filename = "gsva_StromalScore_RawC6.png", plot = p, dpi = 600, width = 3, height = 3)

###############################################################################
#'                             Manuscript: Figure 5g                   '#
###############################################################################                           

########################################################################
gene_data_Fib <- read.csv("~/scRNAseqData of HCC/Fib_Markers.csv") 
filtered_gene_data <- gene_data_Fib %>%
  filter(`Adjusted.P.value` < 0.01) %>%
  arrange(desc(`Log2.Fold.Change`))
top_genes_by_cluster <- filtered_gene_data %>%
  group_by(Cluster) %>%
  slice_head(n = 50) %>%
  summarise(Genes = list(Gene), .groups = 'drop')
gene_sets_list_Fib <- setNames(top_genes_by_cluster$Genes, top_genes_by_cluster$Cluster)
head(gene_sets_list_Fib)

# Creating a new list 'Marker_list' from the specified subsets
Marker_list <- list(
  CAF = gene_sets_list_Fib$CAF,
  Neu_C_1 = DscMarker_list$Neu_C_1,
  Raw_C_2 = RawMarker_list$Raw_C_2,
  Raw_C_6 = RawMarker_list$Raw_C_6
)
########################################################################
# To view the new list
print(Marker_list)
P1T <- readRDS('P1T.rds')

cluster_colors <- c(
  "HMGB2 malignant hepatocyte" = "#7a6cd2",  
  "malignant hepatocyte" = "#70d3ff",   
  "SPP1_Macrophage/CAF" = "#c62d51",   
  "MAP3K12 malignant hepatocyte" = "#60539e",   
  "hepatocyte" = "#c4eba5",   
  "Myofibroblast_Pericyte" = "#d985b4",   
  "Immune_Fibroblast" = "#f1dc92",   
  "Fibroblast" = "#f4f8f9",
  "Myeloid" = "#f1945a"
)


plot <- SpatialDimPlot(P1T,
                       group.by = "DefineTypes",
                       pt.size.factor = 1.5,
                       image.alpha=0.5,
                       cols = cluster_colors, 
                       #stroke =NA, 
                       images = "H&E") +
  theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 16)))

print(plot)
ggsave(filename = "SpatialDimP1T.png", plot = plot, dpi = 600, width = 5, height = 5)


subfolder_path <- paste0('~/Fig5/',"GSVA",'/P1T/RAWEXP/')

if (!dir.exists(subfolder_path)) {
  dir.create(subfolder_path, recursive = TRUE)
}

table(P1T$Dsc_cluster)


########################################################################
expr_matrix <- GetAssayData(P1T, assay = "RNA", slot = "data")
gsva_scores <- gsva(expr_matrix, Marker_list, kcdf="Gaussian",method = "gsva", verbose = FALSE, parallel.sz = 23)
gsva_scores_df <- as.data.frame(t(gsva_scores))
rownames(gsva_scores_df) <- colnames(expr_matrix)
colnames(gsva_scores_df)
P1T <- AddMetaData(P1T, metadata = gsva_scores_df)

for (feature in rownames(gsva_scores)) {
  plot <- SpatialPlot(P1T,
                      features = feature,
                      image.alpha=0.5,
                      images = "H&E",
                      stroke =NA,
                      pt.size.factor = 1.4, alpha = c(0.5, 1.0)) +
    scale_fill_gradientn(colors = c("#f4f6f5",'#fbdfcc',"#eead9b", "#ca4544",'#b91b2f','#842525')) + 
    hide_axis +    
    theme(
      legend.title = element_text(size = 18), 
      text = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 7),
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank()   
    )
  GSVA <- paste0(subfolder_path,'RawExp_', feature, ".png")
  
  ggsave(GSVA, plot, width = 15, height = 12, units = "cm")
}

##############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
########################################################################
subfolder_path <- paste0('~/Fig5/',"GSVA",'/P1T/DSCEXP/')

if (!dir.exists(subfolder_path)) {
  dir.create(subfolder_path, recursive = TRUE)
}

expr_matrix <- GetAssayData(P1T, assay = "DSC", slot = "data")
gsva_scores <- gsva(expr_matrix, Marker_list, kcdf="Gaussian",method = "gsva", verbose = FALSE, parallel.sz = 16)
gsva_scores_df <- as.data.frame(t(gsva_scores))
rownames(gsva_scores_df) <- colnames(expr_matrix)
colnames(gsva_scores_df)
P1T <- AddMetaData(P1T, metadata = gsva_scores_df)

for (feature in rownames(gsva_scores)) {
  plot <- SpatialPlot(P1T,
                      features = feature,
                      image.alpha=0.5,
                      images = "H&E",
                      stroke =NA,
                      pt.size.factor = 1.4, alpha = c(0.5, 1.0)) +
    scale_fill_gradientn(colors = c("#f4f6f5",'#fbdfcc',"#eead9b", "#ca4544",'#b91b2f','#842525'),limits = c(0, 0.6),
                         oob = squish) + 
    hide_axis +    
    theme(
      legend.title = element_text(size = 18), #
      text = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 7),
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank()
    )

  GSVA <- paste0(subfolder_path,'DscM_', feature, ".png")

  ggsave(GSVA, plot, width = 15, height = 12, units = "cm")
}


###############################################################################
#'                             Manuscript: Figure 5h                   '#
###############################################################################

P7T <- readRDS('P7T.rds')

cluster_colors <- c(
  "Malignant_hepatocyte/Myeloid" = "#7a6cd2",  
  "Malignant hepatocyte" = "#70d3ff",   
  "SPP1_Macrophage/CAF" = "#c62d51",   
  "Malignant_hepatocyte/Malignant_cholangiocyte" = "#60539e",   
  "Hepatocyte" = "#c4eba5",   
  "Myofibroblast_Pericyte" = "#d985b4",   
  "Fibroblst_Immune" = "#f1dc92",   
  "Fibroblast" = "#f4f8f9",
  "Myeloid" = "#f1945a",
  "Malignnat_hepatocyte/Fibroblast" = '#3f7fb5'
)


plot <- SpatialDimPlot(P7T,
                       group.by = "DefineTypes",
                       pt.size.factor = 1.5,
                       image.alpha=0.5,
                       cols = cluster_colors, 
                       #stroke =NA, 
                       images = "H&E") +
  theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 16)))

print(plot)
ggsave(filename = "SpatialDimP7T.png", plot = plot, dpi = 600, width = 5, height = 5)

########################################################################
########################################################################

subfolder_path <- paste0('~/Fig5/',"GSVA",'/P7T/RAWEXP/')
if (!dir.exists(subfolder_path)) {
  dir.create(subfolder_path, recursive = TRUE)
}

expr_matrix <- GetAssayData(P7T, assay = "RNA", slot = "data")
gsva_scores <- gsva(expr_matrix, Marker_list, kcdf="Gaussian",method = "gsva", verbose = FALSE, parallel.sz = 16)
gsva_scores_df <- as.data.frame(t(gsva_scores))
rownames(gsva_scores_df) <- colnames(expr_matrix)
colnames(gsva_scores_df)
P7T <- AddMetaData(P7T, metadata = gsva_scores_df)

for (feature in rownames(gsva_scores)) {
  plot <- SpatialPlot(P7T,
                      features = feature,
                      image.alpha=0.5,
                      images = "H&E",
                      stroke =NA,
                      pt.size.factor = 1.4, alpha = c(0.5, 1.0)) +
    scale_fill_gradientn(colors = c("#f4f6f5",'#fbdfcc',"#eead9b", "#ca4544",'#b91b2f','#842525')) + 
    hide_axis +    
    theme(
      legend.title = element_text(size = 18),
      text = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 7),
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank()   
    )
  GSVA <- paste0(subfolder_path,'DscM_', feature, ".png")
  
  ggsave(GSVA, plot, width = 15, height = 12, units = "cm")
}



#############################################################################################
#############################################################################################
#############################################################################################
########################################################################
subfolder_path <- paste0('~/Fig5/',"GSVA",'/P7T/DSCEXP/')
if (!dir.exists(subfolder_path)) {
  dir.create(subfolder_path, recursive = TRUE)
}

expr_matrix <- GetAssayData(P7T, assay = "DSC", slot = "data")
gsva_scores <- gsva(expr_matrix, Marker_list, kcdf="Gaussian",method = "gsva", verbose = FALSE, parallel.sz = 16)
gsva_scores_df <- as.data.frame(t(gsva_scores))
rownames(gsva_scores_df) <- colnames(expr_matrix)
colnames(gsva_scores_df)
P7T <- AddMetaData(P7T, metadata = gsva_scores_df)

if (!dir.exists(subfolder_path)) {
  dir.create(subfolder_path, recursive = TRUE)
}

for (feature in rownames(gsva_scores)) {
  plot <- SpatialPlot(P7T,
                      features = feature,
                      image.alpha=0.5,
                      images = "H&E",
                      stroke =NA,
                      pt.size.factor = 1.4, alpha = c(0.5, 1.0)) +
    scale_fill_gradientn(colors = c("#f4f6f5",'#fbdfcc',"#eead9b", "#ca4544",'#b91b2f','#842525')) + 
    hide_axis +    
    theme(
      legend.title = element_text(size = 18), 
      text = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 7),
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank()  
    )
  GSVA <- paste0(subfolder_path,'DscM_', feature, ".png")
  
  ggsave(GSVA, plot, width = 15, height = 12, units = "cm")
}


###############################################################################
#'                             Manuscript: Figure 5I   #P1T                '#
############################################################################### 

cluster_colors <- c(
  "HMGB2 malignant hepatocyte" = "#7a6cd2",  
  "malignant hepatocyte" = "#70d3ff",   
  "SPP1_Macrophage/CAF" = "#c62d51",   
  "MAP3K12 malignant hepatocyte" = "#60539e",   
  "hepatocyte" = "#c4eba5",   
  "Myofibroblast_Pericyte" = "#d985b4",   
  "Immune_Fibroblast" = "#f1dc92",   
  "Fibroblast" = "#f4f8f9",
  "Myeloid" = "#f1945a"
)

data_to_plot <- FetchData(P1T, vars = c("DefineTypes","Neu_C_1"))

data_long <- reshape2::melt(data_to_plot, id.vars = "DefineTypes")

plot <- ggplot(data_long, aes(x = DefineTypes, y = value, fill = DefineTypes)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  facet_wrap(~ variable, scales = "free_y") +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  scale_fill_manual(values = cluster_colors) +
  labs(title = "", x = "DefineTypes", y = "Expression Score") +
  theme_classic() +
  theme(text = element_text(color = "black", face = "bold"), 
        axis.title = element_text(size = 12,color = "black", face = "bold"),  
        axis.text = element_text(size = 10,color = "black", face = "bold"), 
        strip.text = element_text(size = 10,color = "black", face = "bold"),
        axis.text.x=element_text(angle = 60, hjust = 1,vjust=1),
        legend.position="right"
  ) 

print(plot)
ggsave(filename = "fig51NeuC1.png", plot = plot, dpi = 600, width = 10, height = 6)

########################################## CAF #####################################
data_to_plot <- FetchData(P1T, vars = c("DefineTypes","CAF"))

data_long <- reshape2::melt(data_to_plot, id.vars = "DefineTypes")

plot <- ggplot(data_long, aes(x = DefineTypes, y = value, fill = DefineTypes)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  facet_wrap(~ variable, scales = "free_y") +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  scale_fill_manual(values = cluster_colors) +
  labs(title = "", x = "DefineTypes", y = "Expression Score") +
  theme_classic() +
  theme(text = element_text(color = "black", face = "bold"),  
        axis.title = element_text(size = 12,color = "black", face = "bold"),  
        axis.text = element_text(size = 10,color = "black", face = "bold"),  
        strip.text = element_text(size = 10,color = "black", face = "bold"),
        axis.text.x=element_text(angle = 60, hjust = 1,vjust=1),
        legend.position="right"
  )  

print(plot)
ggsave(filename = "fig51caf.png", plot = plot, dpi = 600, width = 10, height = 6)

###############################################################################
#'                             Manuscript: Figure 5j                '#
############################################################################### 

################# P1T ####################
data_to_plot <- FetchData(P1T, vars = c("DefineTypes", "Neu_C_1", "CAF"))
df1 = as.data.frame(cbind(data_to_plot$Neu_C_1, data_to_plot$CAF))
names(df1) <- c("Neu_C_1", "CAF")

corT = cor.test(df1$Neu_C_1, df1$CAF, method = "spearman")
cor = corT$estimate
pValue = corT$p.value

p1 = ggplot(df1, aes(Neu_C_1, CAF)) +
  geom_point(color = "#ebced2", size = 2) +  
  geom_smooth(method = "lm", formula = y ~ x, color = "#1f5487") +  
  theme_bw() +
  stat_cor(method = 'spearman', aes(label = paste0("Spearman Correlation: ", round(cor, 2), "\nP-value: ", format.pval(pValue, digits = 2))))


p2 = ggMarginal(p1, type = "density", xparams = list(fill = "#e396bb"), yparams = list(fill = "#5ba2c2"))
print(p2)
ggsave("Fig5jp1t.png", plot = p2, width = 4.2, height = 4, dpi = 600)


################# P7T ####################
data_to_plot <- FetchData(P7T, vars = c("DefineTypes", "Neu_C_1", "CAF"))
df1 = as.data.frame(cbind(data_to_plot$Neu_C_1, data_to_plot$CAF))
names(df1) <- c("Neu_C_1", "CAF")

corT = cor.test(df1$Neu_C_1, df1$CAF, method = "spearman")
cor = corT$estimate
pValue = corT$p.value

p1 = ggplot(df1, aes(Neu_C_1, CAF)) +
  geom_point(color = "#ebced2", size = 2) +  
  geom_smooth(method = "lm", formula = y ~ x, color = "#1f5487") + 
  theme_bw() +
  stat_cor(method = 'spearman', aes(label = paste0("Spearman Correlation: ", round(cor, 2), "\nP-value: ", format.pval(pValue, digits = 2))))


p2 = ggMarginal(p1, type = "density", xparams = list(fill = "#e396bb"), yparams = list(fill = "#5ba2c2"))
print(p2)
ggsave("Fig5jp7t.png", plot = p2, width = 4.2, height = 4, dpi = 600)
