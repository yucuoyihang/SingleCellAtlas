### Figure 4G

library(Seurat)

cols_celltype = c('CD4T' = '#8DD3C7',
                  'CD8T' = '#FFFFB3',
                  'B' = '#FB8072',
                  'NK' = '#BEBADA',
                  'DC' = '#B3DE69',
                  'Plasma' = '#80B1D3',
                  'Monocyte' = '#FDB462')

seurat_obj = readRDS("valid/Valid_PBMC_CSF.rds")

pdf(file = "Figure 4G_umap.pdf", width = 5.5, height = 4)
DimPlot(seurat_obj, group.by="celltype", cols = cols_celltype, raster = T, split.by = "origin")+NoAxes()
dev.off()

cols_celltype_monocyte = c('cMono' = '#FFB7B7',
                  'ncMono_LYVE1' = '#ADE892',
                  'ncMono' = '#B5CEEA')

seurat_obj = readRDS("valid/Valid_PBMC_CSF.rds")

pdf(file = "Figure 4G_umap_monocyte.pdf", width = 4, height = 4)
DimPlot(seurat_obj, group.by="fine", cols = cols_celltype_monocyte, raster = T, split.by = "origin")+NoAxes()
dev.off()


### Figure 4H

library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(patchwork)

palette_volcano <- c("Up" = "#8B0000", "Down" = "#4682B4", "NoSig" = "#B3B3B3")

highlight_genes <- read.csv("DEGs_CSF_Biological_AD vs CN_A-_fine.csv", row.names = 1) %>%
  filter(celltype == "ncMono_05_LYVE1", p_val_adj < 0.05, abs(avg_log2FC) > 1.5) %>%
  pull(gene)

process_and_plot <- function(file_path, plot_title, target_genes) {
  dt <- read.csv(file_path, row.names = 1) %>%
    rownames_to_column("gene") %>%
    filter(!str_detect(gene, "^MT-|^RPL|^RPS|^HLA-")) %>%
    mutate(
      status = case_when(
        avg_log2FC > 0.25 & p_val_adj < 0.05 ~ "Up",
        avg_log2FC < -0.25 & p_val_adj < 0.05 ~ "Down",
        TRUE ~ "NoSig"
      ),
      status = factor(status, levels = c("Up", "Down", "NoSig"))
    )
  
  hl <- dt %>% filter(gene %in% target_genes)
  
  p_volcano <- ggplot(dt, aes(x = avg_log2FC, y = -log10(p_val_adj), color = status)) +
    geom_point(alpha = 0.9, size = 2.2) +
    geom_point(data = hl, shape = 1, stroke = 0.5, size = 2.5, color = "black", alpha = 1) +
    geom_hline(yintercept = 1.3, color = "grey40", linetype = "longdash", linewidth = 0.25) +
    geom_vline(xintercept = c(-0.25, 0.25), color = "grey40", linetype = "longdash", linewidth = 0.25) +
    scale_color_manual(values = palette_volcano, drop = FALSE) +
    labs(x = "Log2FC", y = "-Log10P-adj", title = plot_title) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 10, color = "black")
    )
  
  return(p_volcano)
}

plot <- process_and_plot("Valid_AD_DEGs_ncMono_LYVE1.csv", "ncMono_LYVE1  AD−DEGs", highlight_genes)

ggsave("Figure 4H.pdf", plot = plot, width = 3, height = 2)


### Figure 4J

library(Seurat)
library(dplyr)
library(clusterProfiler)
library(enrichplot)

seurat_obj = readRDS("valid/Human_Myeloid_atlas_v2.RDS")

valid_regions <- c("Spinal_cord")
valid_diseases <- c("AD_unspecified_Disease", "AD_early_Disease", "AD_late_Disease", "Control_Control")
valid_subclasses <- c("Microglia", "CAM", "MDM")

seurat_obj <- subset(
  seurat_obj,
  subset = !(Region %in% valid_regions) & 
    Disease.condition %in% valid_diseases & 
    Subclass_label_final %in% valid_subclasses
)

Idents(seurat_obj) <- "Subclass_label_final"

target_gene_set <- read.csv("AD_ncMono_05_LYVE1_marker.csv", row.names = 1) %>%
  rownames_to_column("gene") %>%
  filter(p_val_adj < 0.05, avg_log2FC > 2) %>%
  pull(gene)

term2gene <- tibble(
  term = "Target_Subpop_Markers",
  gene = target_gene_set
)

celltypes <- c("Microglia", "CAM", "MDM")

gsea_results_list <- map(celltypes, function(ct) {
  deg <- FindMarkers(
    object = seurat_obj,
    ident.1 = ct,
    min.pct = 0.1,
    logfc.threshold = 0,
    verbose = FALSE
  )
  
  ranked_list <- deg %>%
    rownames_to_column("gene") %>%
    arrange(desc(avg_log2FC)) %>%
    deframe()
  
  GSEA(
    geneList = ranked_list,
    TERM2GENE = term2gene,
    pvalueCutoff = 1,
    minGSSize = 5,
    maxGSSize = 5000,
    eps = 0,
    verbose = FALSE
  )
}) %>% 
  set_names(celltypes)

summary_table <- map_dfr(gsea_results_list, ~ as.data.frame(.x), .id = "CellType_B")

summary_table %>%
  select(CellType_B, NES, pvalue, p.adjust) %>%
  print()

p_gsea <- gseaplot2(
  gsea_results_list[["CAM"]],
  geneSetID = 1,
  color = "#8A1619",
  rel_heights = c(0.5, 0.2, 0.1)
)

ggsave("Figure_4J.pdf", plot = p_gsea, width = 4, height = 3)





### Figure 4K

library(Seurat)
library(dplyr)
library(clusterProfiler)
library(enrichplot)

seurat_obj <- readRDS("valid/Human_Myeloid_atlas_v2.RDS")

valid_diseases <- c("AD_unspecified_Disease", "AD_early_Disease", "AD_late_Disease", "Control_Control")
valid_subclasses <- c("Microglia", "CAM", "MDM")

seurat_obj <- subset(
  seurat_obj,
  subset = !(Region %in% "Spinal_cord") & 
    Disease.condition %in% valid_diseases & 
    Subclass_label_final %in% valid_subclasses
)

seurat_obj$Disease.condition <- factor(
  seurat_obj$Disease.condition,
  levels = c("Control_Control", "AD_early_Disease", "AD_late_Disease", "AD_unspecified_Disease"),
  labels = c("CN", "AD", "AD", "AD")
)

target_genes <- c("C1QC", "C1QB", "C1QA", "LYVE1", "C3", "APOE", "ANKAR", "CSF1R", "AXL")
seurat_obj <- AddModuleScore(object = seurat_obj,  features = list(target_genes),  ctrl = 100,  name = "gs")
gradient_cols <- c("lightblue3", "lightblue", "white", "red", "red4")


p1 <- DotPlot(
  object = seurat_obj,
  features = target_genes,
  group.by = "Disease.condition",
  idents = "CAM"
) +
  scale_color_gradientn(colours = gradient_cols) +
  RotatedAxis() +
  labs(x = NULL, y = NULL, title = "CAM Marker Genes") +
  theme(axis.text = element_text(color = "black"))

p2 <- DotPlot(
  object = seurat_obj,
  features = "gs1",
  group.by = "Disease.condition",
  idents = "CAM"
) +
  scale_color_gradientn(colours = gradient_cols) +
  RotatedAxis() +
  labs(x = NULL, y = NULL, title = "Score") +
  theme(axis.text = element_text(color = "black"))

p2$data$pct.exp <- 100
p2 <- p2 + 
  scale_size(range = c(6, 6)) + 
  guides(size = "none")

final_plot <- p1 + p2 + plot_layout(ncol = 2, widths = c(1, 0.2))

ggsave("Figure 4K.pdf", plot = final_plot, width = 10, height = 2.5)




### Figure 4N

library(Seurat)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(tidyverse)

palette_group <- c("AD" = "#CB3B3B", "CN" = "#5D8BB5")
target_genes  <- c("C1QC", "C1QA", "LYVE1", "C3AR1")

seurat_obj <- readRDS("valid/Immune_cells_with_metadata.rds")

deg_genes <- read.csv("AD_ncMono_05_LYVE1_marker.csv", row.names = 1) %>%
  filter(avg_log2FC > 1, p_val_adj < 0.05) %>%
  rownames()

Idents(seurat_obj) <- "cell_type_high_resolution" 
seurat_obj <- subset(seurat_obj, idents = "CAMs")
seurat_obj <- AddModuleScore(object = seurat_obj, features = list(deg_genes), ctrl = 100, name = "gs")

clean_meta <- seurat_obj@meta.data %>%
  mutate(Braak = as.numeric(as.character(Braak))) %>%
  filter(
    !(group == "AD" & Braak == 1),
    !(group == "CN" & Braak == 4)
  )

p_score <- ggplot(clean_meta, aes(x = Braak, y = gs1)) +
  geom_smooth(aes(color = group, fill = group), method = "lm", se = TRUE, span = 1) +
  stat_cor(aes(color = group), method = "spearman", size = 4, label.x = 0.5, label.y = c(0.03, 0.034)) +
  scale_color_manual(values = palette_group) +
  scale_fill_manual(values = palette_group) +
  scale_x_continuous(breaks = 0:6) +
  labs(x = "Braak Stage", y = NULL, title = "Brain CAMs") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    legend.position = "none"
  )

ggsave("Figure 4N_signature.pdf", plot = p_score, width = 3.75, height = 3.5)

plot_gene_data <- FetchData(seurat_obj, vars = c("Braak", "group", target_genes), cells = rownames(clean_meta)) %>%
  mutate(Braak = as.numeric(as.character(Braak))) %>%
  pivot_longer(cols = all_of(target_genes), names_to = "Gene", values_to = "Expression")

p_genes <- ggplot(plot_gene_data, aes(x = Braak, y = Expression)) +
  geom_smooth(aes(color = group, fill = group), method = "loess", se = TRUE, span = 1) +
  facet_wrap(~ Gene, nrow = 2, scales = "free_y") +
  scale_color_manual(values = palette_group) +
  scale_fill_manual(values = palette_group) +
  scale_x_continuous(breaks = 0:6) +
  labs(x = "Braak Stage", y = "Expression", title = "Brain CAMs") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )

ggsave("Figure 4N_genes.pdf", plot = p_genes, width = 5, height = 5)


