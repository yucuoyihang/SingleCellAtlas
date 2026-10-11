### Figure 6B

library(Seurat)
library(tidyverse)
library(edgeR)

celltype_colors <- c(
  'CD4T_Tn_01_GAS5'='#e5d25b', 'CD4T_Tn_02_SELL'='#599014', 'CD4T_Tn_03_LEF1'='#cd5d1a', 
  'CD4T_Tn_04_TCF7'='#a82d06', 'CD4T_Tn_05_IL7R'='#4592BF', 'CD4T_Tcm_01_CD69'='#d38219', 
  'CD4T_Tcm_02_GPR183'='#74a764', 'CD4T_Tcm_03_TCF7'='#8ca2b4', 'CD4T_Tcm_04_S100A4'='#cbb190', 
  'CD4T_Tem'='#e7ca8d', 'CD4T_CTL'='#b03e69', 'CD4T_Treg'='#9d9ec3', 'CD8T_Tn_01_IL7R'='#593202', 
  'CD8T_Tn_02_TCF7'='#467599', 'CD8T_Tn_03_LEF1'='#d36135', 'CD8T_Tcm'='#ebafa4', 
  'CD8T_Tem_01_CX3CR1'='#5e8a89', 'CD8T_Tem_02_GZMK'='#ecd577', 'CD8T_Teff_01_GZMH'='#7c606c', 
  'CD8T_Teff_02_KLRG1'='#5c6489', 'CD8T_Tex'='#602a2a', 'NK_CD56_01_CCL4'='#6fadd3', 
  'NK_CD56_02_XCL1'='#1c909a', 'NK_CD16_01_CD44'='#fe7e02', 'NK_CD16_02_PRF1'='#e4191c', 
  'NK_CD16_03_IGFBP7'='#397eb7', 'NK_CD16_04_KLRC2'='#4D377A', 'NK_CD16_05_CX3CR1'='#E5E7F3', 
  'NK_CD16_06_CMC1'='#68C3A6', 'Bn_01_TCL1A'='#6a7fe1', 'Bn_02_IL4R'='#c5aecb', 
  'Bn_03_SAMHD1'='#ca7474', 'Bm_01_SLC25A6'='#F1EFDE', 'Bm_02_SAMD3'='#8A1619', 
  'Bm_03_IFITM3'='#00589C', 'Batyp_01_ITPR1'='#d18e95', 'Batyp_02_FCRLA'='#b5b47b', 
  'Bact'='#edf2bc', 'Bgc'='#FFB7BE', 'Bcyc'='#87C48F', 'Plasma'='#B3DE69', 
  'cMono_01_LRMDA'='#234862', 'cMono_02_S100A12'='#9078d4', 'cMono_03_TREM1'='#A71D1D', 
  'cMono_04_JUN'='#342488', 'ncMono_01_CALHM6'='#797D7D', 'ncMono_02_TBC1D8'='#c5cc95', 
  'ncMono_03_CXCR4'='#40032b', 'ncMono_04_BCL2'='#dbb679', 'ncMono_05_LYVE1'='#c43c40', 
  'cycMono'='#4E4A8C', 'cDC1'='#115F5F', 'cDC2'='#238b5c', 'cDC2_stress'='#59b5b2', 
  'pDC'='#9fcccf', 'tDC'='#008BBD', 'MicDC'='#FF9F80', 'Cycling'='#FCCDE5'
)

seurat_obj <- readRDS("AD_all_paired.rds")

seurat_obj_csf <- subset(seurat_obj, subset = origin == "CSF")

meta_csf <- seurat_obj_csf@meta.data %>%
  group_by(orig.ident, celltype) %>%
  mutate(celltype_total = n()) %>%
  group_by(orig.ident, celltype, fine) %>%
  mutate(fine_count = n()) %>%
  ungroup() %>%
  mutate(ratio = fine_count / celltype_total) %>%
  filter(fine == "ncMono_05_LYVE1") %>%
  distinct(orig.ident, ratio) %>%
  mutate(patient = str_split_i(orig.ident, "_", 1)) %>%
  rename(csf_ncMono_05_LYVE1_ratio = ratio)

seurat_obj_pbmc <- subset(seurat_obj, subset = origin == "Blood")
pbmc_celltypes <- unique(seurat_obj_pbmc$fine)

corr_results_list <- map(pbmc_celltypes, function(ct) {
  seurat_ct <- subset(seurat_obj_pbmc, subset = fine == ct)
  
  agg_counts <- AggregateExpression(
    seurat_ct, 
    group.by = "orig.ident", 
    assays = "RNA", 
    layers = "counts"
  )$RNA
  
  agg_counts_mat <- as.matrix(agg_counts)
  
  if (ncol(agg_counts_mat) < 5) {
    message("Low sample size and pass: ", ct)
    return(NULL)
  }
  
  log2_cpm <- cpm(agg_counts_mat, log = TRUE, prior.count = 1)
  
  avg_exp_long <- as.data.frame(log2_cpm) %>%
    rownames_to_column("gene") %>%
    pivot_longer(cols = -gene, names_to = "orig.ident", values_to = "expression") %>%
    mutate(patient = str_split_i(orig.ident, "_", 1))
  
  merged_data <- inner_join(avg_exp_long, meta_csf, by = "patient")
  paired_sample_count <- n_distinct(merged_data$patient)
  
  merged_data %>%
    group_by(gene) %>%
    filter(sd(expression) > 0) %>%
    summarise(
      cor_coef = cor(expression, csf_ncMono_05_LYVE1_ratio, method = "spearman"),
      p_value = cor.test(expression, csf_ncMono_05_LYVE1_ratio, method = "spearman")$p.value,
      sample_size = paired_sample_count,
      .groups = "drop"
    ) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value)
}) %>% 
  set_names(pbmc_celltypes) %>% 
  compact() 

final_corr_df <- bind_rows(corr_results_list, .id = "celltype")

plot_data <- final_corr_df %>%
  filter(p_value < 0.05, cor_coef > 0.5) %>%
  count(celltype, name = "Freq") %>%
  filter(Freq > 5) %>%
  mutate(celltype = reorder(celltype, -Freq))

p_bar <- ggplot(plot_data, aes(x = celltype, y = Freq, fill = celltype)) +
  geom_col() +
  scale_fill_manual(values = celltype_colors) +
  scale_y_sqrt() +
  labs(
    x = "Cell types", 
    y = "Number of ncMono_05_LYVE1-positively correlated genes"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    legend.position = "none"
  )

ggsave("Figure 6B.pdf", plot = p_bar, width = 6, height = 5)


### Figure 6C

library(nichenetr)
library(Seurat)
library(tidyverse)
library(ggplot2)

lr_network <- readRDS("resources/nichenet/lr_network_human_21122021.rds")
ligand_target_matrix <- readRDS("resources/nichenet/ligand_target_matrix_nsga2r_final.rds")
weighted_networks <- readRDS("resources/nichenet/weighted_networks_nsga2r_final.rds")

seurat_obj <- readRDS("AD_all_paired.rds")
seurat_obj@meta.data$label = ifelse(seurat_obj@meta.data$fine == "ncMono_05_LYVE1", "ncMono_05_LYVE1", seurat_obj@meta.data$celltype)
Idents(seurat_obj) = as.factor(seurat_obj$fine)

receiver = c("ncMono_05_LYVE1")
expressed_genes_receiver = get_expressed_genes(receiver, seurat_obj, pct = 0.05)

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

sender_celltypes <- c("NK_CD16_05_CX3CR1")
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_obj, 0.01)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

condition_oi <-  "AD"
condition_reference <- "CN"
seurat_obj_receiver <- subset(seurat_obj, idents = receiver)
DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, group.by = "group", min.pct = 0.01) %>% rownames_to_column("gene")
geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
ligand_activities <- predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands_focused)
ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities
best_upstream_ligands <- ligand_activities %>% top_n(20, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)

vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>% column_to_rownames("test_ligand")%>% select(aupr)%>% arrange(aupr) %>% as.matrix(ncol = 1)
vis_ligand_aupr = as.data.frame(vis_ligand_aupr)
min_val <- min(vis_ligand_aupr$aupr)
max_val <- max(vis_ligand_aupr$aupr)
mid_val <- 0.29 
mid_prop <- (mid_val - min_val) / (max_val - min_val)
p1=make_heatmap_ggplot(as.matrix(vis_ligand_aupr), "Prioritized ligands", "Ligand activity", legend_title = "AUPR", color = "darkorange") + theme(axis.text.x.top = element_blank()) + scale_fill_gradientn(colors = c("#F7F7F7","#92A9D8", "#2670BB"), values = c(0, mid_prop, 1))

ligand_receptor_links_df <- get_weighted_ligand_receptor_links(best_upstream_ligands, expressed_receptors, lr_network, weighted_networks$lr_sig) 
vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(ligand_receptor_links_df, best_upstream_ligands, order_hclust = "both") 
p2=(make_heatmap_ggplot(t(vis_ligand_receptor_network), y_name = "Ligands", x_name = "Receptors", color = "#5C3A93", legend_title = "Prior interaction potential"))

pdf(file=paste0("Figure 6C.pdf"),width=7, height=4.5)
p1+p2+plot_layout(width=c(1,15))
dev.off()

