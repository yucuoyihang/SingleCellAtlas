### Figure 3A

library(Seurat)

seurat_obj = readRDS("AD_Mono.rds")

pdf(file = "Figure 3A_umap.pdf", width = 5.5, height = 3.75)
DimPlot(seurat_obj, group.by="fine", raster=F, cols = c("#234862", "#9078d4", "#A71D1D", "#342488", "#4E4A8C", "#dbb679", "#c43c40", "#797D7D", "#c5cc95", "#40032b"))+NoAxes()
dev.off()

library(ggplot2)
library(dplyr)

cols_fine = c('CD4T_Tn_01_GAS5' = '#e5d25b',
              'CD4T_Tn_02_SELL' = '#599014',
              'CD4T_Tn_03_LEF1' = '#cd5d1a',
              'CD4T_Tn_04_TCF7' = '#a82d06',
              'CD4T_Tn_05_IL7R' = '#4592BF',
              'CD4T_Tcm_01_CD69' = '#d38219',
              'CD4T_Tcm_02_GPR183' = '#74a764',
              'CD4T_Tcm_03_TCF7' = '#8ca2b4',
              'CD4T_Tcm_04_S100A4' = '#cbb190',
              'CD4T_Tem' = '#e7ca8d',
              'CD4T_CTL' = '#b03e69',
              'CD4T_Treg' = '#9d9ec3',
              'CD8T_Tn_01_IL7R' = '#593202',
              'CD8T_Tn_02_TCF7' = '#467599',
              'CD8T_Tn_03_LEF1' = '#d36135',
              'CD8T_Tcm' = '#ebafa4',
              'CD8T_Tem_01_CX3CR1' = '#5e8a89',
              'CD8T_Tem_02_GZMK' = '#ecd577',
              'CD8T_Teff_01_GZMH' = '#7c606c',
              'CD8T_Teff_02_KLRG1' = '#5c6489',
              'CD8T_Tex' = '#602a2a',
              'NK_CD56_01_CCL4' = '#6fadd3',
              'NK_CD56_02_XCL1' = '#1c909a',
              'NK_CD16_01_CD44' = '#fe7e02',
              'NK_CD16_02_PRF1' = '#e4191c',
              'NK_CD16_03_IGFBP7' = '#397eb7',
              'NK_CD16_04_KLRC2' = '#4D377A',
              'NK_CD16_05_CX3CR1' = '#E5E7F3',
              'NK_CD16_06_CMC1' = '#68C3A6',
              'Bn_01_TCL1A' = '#6a7fe1',
              'Bn_02_IL4R' = '#c5aecb',
              'Bn_03_SAMHD1' = '#ca7474',
              'Bm_01_SLC25A6' = '#F1EFDE',
              'Bm_02_SAMD3' = '#8A1619',
              'Bm_03_IFITM3' = '#00589C',
              'Batyp_01_ITPR1' = '#d18e95',
              'Batyp_02_FCRLA' = '#b5b47b',
              'Bact' = '#edf2bc',
              'Bgc' = '#FFB7BE',
              'Bcyc' = '#87C48F',
              'Plasma' = '#B3DE69',
              'cMono_01_LRMDA' = '#234862',
              'cMono_02_S100A12' = '#9078d4',
              'cMono_03_TREM1' = '#A71D1D',
              'cMono_04_JUN' = '#342488',
              'ncMono_01_CALHM6' = '#797D7D',
              'ncMono_02_TBC1D8' = '#c5cc95',
              'ncMono_03_CXCR4' = '#40032b',
              'ncMono_04_BCL2' = '#dbb679',
              'ncMono_05_LYVE1' = '#c43c40',
              'cycMono' = '#4E4A8C',
              'cDC1' = '#115F5F',
              'cDC2' = '#238b5c',
              'cDC2_stress' = '#59b5b2',
              'pDC' = '#9fcccf',
              'tDC' = '#008BBD',
              'MicDC' = '#FF9F80',
              'Cycling' = '#FCCDE5')

deg_counts <- read.csv("DEGs_CSF_Biological_AD vs CN_A-_fine.csv") %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5) %>%
  mutate(direc = if_else(avg_log2FC > 0, "up", "down")) %>%
  count(celltype, direc, name = "Freq") %>%
  group_by(celltype) %>%
  mutate(Total = sum(Freq)) %>%
  ungroup() %>%
  mutate(
    celltype = reorder(celltype, -Total),
    LogSize = log10(Freq + 1),
    PlotFreq = if_else(direc == "down", -LogSize, LogSize)
  )

pdf(file = "Figure 3A_bar.pdf", width = 6, height = 3)
ggplot(deg_counts, aes(x = celltype, y = PlotFreq, fill = celltype)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  scale_y_continuous(
    breaks = c(-3, -2, -1, 0, 1, 2, 3),
    labels = c("1000", "100", "10", "0", "10", "100", "1000")
  ) +
  scale_fill_manual(values = cols_fine) +
  labs(x = NULL, y = "Number of DEGs", fill = "Direction") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )
dev.off()




### Figure 3B

library(Seurat)
library(ggplot2)

seurat_obj = readRDS("AD_Mono.rds")

pdf(file = "Figure 3B.pdf", width = 6.5, height = 3.5)
DotPlot(seurat_obj, features = c("C1QC", "C1QB", "C1QA", "A2M", "APOE", "TREM2", "C3", "MKI67"), group.by="fine")+RotatedAxis()+xlab("") + ylab("")+ scale_color_gradientn(colours = c("lightblue3", "lightblue", "white", "red", "red4"))
dev.off()


### Figure 3C

library(Seurat)
library(ggplot2)
library(clusterProfiler)
library(ggpubr)
library(DOSE)
library(stringr)

seurat_obj = readRDS("AD_Mono.rds")

marker = FindMarkers(seurat_obj, ident.1 = c("ncMono_05_LYVE1"), ident.2 = c("cMono_01_LRMDA", "cMono_02_S100A12", "cMono_03_TREM1", "cMono_04_JUN", "ncMono_01_CALHM6", "ncMono_02_TBC1D8", "ncMono_03_CXCR4", "cycMono"), logfc.threshold = 0.5, only.pos = T)
marker = subset(marker, marker$p_val_adj < 0.05)
marker = subset(marker, marker$avg_log2FC > 1)
gene = rownames(marker)
gene <- gene[-grep(pattern = "^MT-|^RPL|^RPS|^HLA-", x = gene)]
gene<-bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
GO<-enrichGO(gene=gene$ENTREZID,OrgDb='org.Hs.eg.db',keyType="ENTREZID",ont="ALL",pAdjustMethod="BH",pvalueCutoff=0.5,qvalueCutoff=0.5,minGSSize=10,maxGSSize=500,readable=TRUE)
write.csv(summary(GO), paste0("GO_ncMono_05_LYVE1.csv"),row.names =F) 


### Figure 3D

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)

seurat_obj = readRDS("AD_Mono.rds")

cell_ratio <- seurat_obj@meta.data %>%
  filter(origin == "CSF") %>%
  mutate(
    group = case_when(
      b_diagnosis %in% c("MCI_due_to_non-AD", "Dementia_due_to_non-AD") ~ "non-AD",
      TRUE ~ b_diagnosis
    ),
    group = factor(group, levels = c("non-AD", "Dementia_due_to_AD", "MCI_due_to_AD", "Preclinical_AD", "CN_A-"))
  ) %>%
  count(fine, group, name = "count") %>%
  group_by(group) %>%
  mutate(Freq = (count / sum(count)) * 100) %>%
  ungroup()

palette_cols <- c("#234862", "#9078d4", "#A71D1D", "#342488", "#4E4A8C", 
                  "#dbb679", "#c43c40", "#797D7D", "#c5cc95", "#40032b")

p1 <- ggplot(cell_ratio, aes(x = Freq, y = group, fill = fine)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = palette_cols) +
  labs(x = "Proportion", y = NULL, title = "Mono in CSF") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "none"
  )

ggsave("Figure 3D.pdf", plot = p1, width = 4.5, height = 5)


### Figure 3E

library(Seurat)
library(ggplot2)

seurat_obj = readRDS("AD_Mono.rds")
seurat_obj = subset(seurat_obj, idents = "ncMono_05_LYVE1")

genes_of_interest = c("OOEP", "APOO", "RBFOX2", "PLEKHA7", "ARMH1", "TEX14", "RNU12", "MT1G", "CMIP", "CHST11", "LRMDA", "CIITA", "LRP1", "LYVE1", "MACRO", "APOE", "GZMK", "TNFAIP3", "TRBC1", "TRAC", "IL7R")

pdf(file = "Figure 3E.pdf", width = 7, height = 3)
DotPlot(seurat_obj, features = genes_of_interest, group.by="b_diagnosis")+RotatedAxis()+xlab("") + ylab("")+ scale_color_gradientn(colours = c("lightblue3", "lightblue", "white", "red", "red4"))
dev.off()



### Figure 3F- GO

library(Seurat)
library(ggplot2)
library(clusterProfiler)
library(ggpubr)
library(DOSE)
library(stringr)

seurat_obj = readRDS("AD_Mono.rds")
seurat_obj = subset(seurat_obj, idents = "ncMono_05_LYVE1")
Idents(seurat_obj) = as.factor(seurat_obj$b_diagnosis)

groups = c("Dementia_due_to_AD", "MCI_due_to_AD", "Preclinical_AD")
for (i in 1:length(groups)){
  marker = FindMarkers(seurat_obj, ident.1 = c(paste0(groups[i])), ident.2 = c("CN_A-"), logfc.threshold = 0.5, min.pct = 0.1, only.pos = T)
  marker = subset(marker, marker$p_val_adj < 0.05)
  marker = subset(marker, marker$avg_log2FC > 0.5)
  gene = rownames(marker)
  gene <- gene[-grep(pattern = "^MT-|^RPL|^RPS|^HLA-", x = gene)]
  gene<-bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  GO<-enrichGO(gene=gene$ENTREZID,OrgDb='org.Hs.eg.db',keyType="ENTREZID",ont="ALL",pAdjustMethod="BH",pvalueCutoff=0.5,qvalueCutoff=0.5,minGSSize=10,maxGSSize=500,readable=TRUE)
  write.csv(summary(GO), paste0("GO_ncMono_05_LYVE1_",groups[i],"_vs_CN.csv"),row.names =F) 
}



### Figure 3G
library(ggplot2)
library(ggpubr)
df <- read.delim(pipe("pbpaste"))
df$T_status <- factor(df$T_status, levels = c("Negative", "Positive"))
my_cols <- c("Negative" = "#4DBBD5", "Positive" = "#E64B35")
p <- ggplot(df, aes(x = T_status, y = ncMono_05_LYVE1)) +
  geom_boxplot(width = 0.4, fill = "white", color = "black", linewidth = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = T_status), width = 0.03, size = 3, alpha = 0.5) +
  scale_color_manual(values = my_cols) +
  labs(title = "P = 0.023", x = "Tau status", y = "CSF ncMono_05_LYVE1 proportion (%)") +
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, color = "black", size = 15), legend.position = "none", axis.text = element_text(color = "black"), axis.title.y = element_text(face = "plain"), axis.line = element_line(linewidth = 0.7))
p
ggsave("ncMono_05_LYVE1_Tau_status.pdf", plot = p, width = 4, height = 5)



### Figure 3H
library(ggplot2)
library(dplyr)
library(ggpubr)
df <- read.table(pipe("pbpaste"), header = TRUE, sep = "\t", check.names = FALSE)
cell_var <- names(df)[1]
v <- "GFAP"
dat <- df[, c(cell_var, v)] 
colnames(dat) <- c("cell", "x")
dat <- dat %>% filter(!is.na(cell), !is.na(x))
med <- median(dat$x, na.rm = TRUE)
dat$group <- factor(ifelse(dat$x >= med, "High", "Low"), levels = c("Low", "High"))
p <- ggplot(dat, aes(x = group, y = cell, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) + geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_y_continuous(limits = c(0, 100)) + theme_bw() +
  labs(title = v, x = paste0("Plasma ", v, " (pg/ml)"), y = "CSF ncMono_05_LYVE1 proportion (%)") +
  scale_fill_manual(values = c("Low" = "#6BAED6", "High" = "#FB6A4A")) + theme(legend.position = "none")
print(p)
ggsave("GFAP_ncMono_05_LYVE1.pdf", plot = p, width = 4, height = 5)
