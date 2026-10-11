# Protein abundance visualization and external validation ======================

library(data.table)
library(dplyr)
library(ggplot2)
library(gghalves)
library(ggpubr)
library(cowplot)

# Colors used for diagnostic groups
group_colors <- c(
  "CN" = "#6090c1", "MCI" = "#f2724d", "Dementia" = "#d6312d", "CI" = "#d6312d",
  "CN A-" = "#6090c1", "Preclinical AD" = "#fde395",
  "MCI due to AD" = "#f2724d", "Dementia due to AD" = "#d6312d"
)

# Half-violin/boxplot with pairwise Wilcoxon tests
plot_half_violin <- function(data, tissue_type, gene, diag_col, unit_label = "") {
  plot_dat <- data %>%
    filter(Tissue == tissue_type, !is.na(.data[[diag_col]]), !is.na(.data[[gene]]))

  groups <- levels(droplevels(plot_dat[[diag_col]]))
  comparisons <- combn(groups, 2, simplify = FALSE)

  sig_comparisons <- Filter(function(comp) {
    wilcox.test(
      as.formula(paste(gene, "~", diag_col)),
      data = plot_dat %>% filter(.data[[diag_col]] %in% comp)
    )$p.value < 0.05
  }, comparisons)

  nudge_val <- -0.15

  ggplot(plot_dat, aes(x = .data[[diag_col]], y = .data[[gene]])) +
    geom_boxplot(
      fill = "white", linewidth = 0.3, outlier.shape = NA, width = 0.5,
      position = position_nudge(x = nudge_val)
    ) +
    geom_half_violin(
      aes(color = .data[[diag_col]], fill = .data[[diag_col]]),
      side = "R", adjust = 1, trim = FALSE, alpha = 0.6, linewidth = 0.4, width = 1
    ) +
    geom_point(
      aes(x = as.numeric(as.factor(.data[[diag_col]])) + nudge_val,
          fill = .data[[diag_col]]),
      shape = 21, stroke = 0.02,
      position = position_jitter(width = 0.05, height = 0, seed = 7),
      alpha = 0.4, size = 1.5
    ) +
    scale_color_manual(values = group_colors) +
    scale_fill_manual(values = group_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.07)),
                       labels = scales::label_comma()) +
    labs(x = "", y = "", title = paste0(tissue_type, " - ", gene, unit_label)) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, color = "black", size = 10),
      axis.text.y = element_text(color = "black", size = 9),
      axis.line = element_line(color = "black", linewidth = 0.35),
      plot.title = element_text(hjust = 0.5, size = 12)
    ) +
    stat_compare_means(
      comparisons = sig_comparisons, method = "wilcox.test",
      label = "p.signif", hide.ns = TRUE, step.increase = 0.1,
      symnum.args = list(
        cutpoints = c(0, 0.01, 0.05, 1),
        symbols = c("**", "*", "ns")
      ),
      tip.length = 0, size = 4.5
    )
}

# Discovery cohort =============================================================

clinical <- fread("HS_protein_raw.tsv") %>%
  mutate(
    Clinical_diagnosis = factor(Clinical_diagnosis, levels = c("CN", "MCI", "Dementia")),
    Biological_diagnosis = factor(
      Biological_diagnosis,
      levels = c("CN_A-", "Preclinical_AD", "MCI_due_to_AD", "Dementia_due_to_AD"),
      labels = c("CN A-", "Preclinical AD", "MCI due to AD", "Dementia due to AD")
    )
  )

p_igfl3_bio <- plot_half_violin(
  clinical, "Blood", "IGFL3", "Biological_diagnosis"
)

p_igfl3_clin <- plot_half_violin(
  clinical, "Blood", "IGFL3", "Clinical_diagnosis"
)

p_spp1_bio <- plot_half_violin(
  clinical, "CSF", "SPP1", "Biological_diagnosis"
)

# ADNI validation: CSF SPP1 ====================================================

dat_ADNI <- fread("CSF_SOMAscan_covariates_CSF708.csv") %>%
  filter(!Diagnosis %in% "Non AD", !is.na(SPP1)) %>%
  mutate(
    Biological_diagnosis = factor(
      Diagnosis,
      levels = c("CN A-", "Preclinical AD", "MCI due to AD", "Dementia due to AD")
    ),
    Tissue = "CSF"
  )

p_spp1_ADNI <- plot_half_violin(
  dat_ADNI, "CSF", "SPP1", "Biological_diagnosis", unit_label = " (ADNI)"
)

# MAP validation: plasma IGFL3 =================================================

dat_MAP <- fread("OhNM2025_ROSMAP_plasma_Soma7k_protein_level_ANML_log10.tsv") %>%
  filter(study == "MAP", Diagnosis %in% c("NCI", "MCI", "AD")) %>%
  mutate(
    Clinical_diagnosis = ifelse(Diagnosis %in% c("MCI", "AD"), "CI", "CN"),
    Clinical_diagnosis = factor(Clinical_diagnosis, levels = c("CN", "CI")),
    Tissue = "Blood"
  )

p_igfl3_MAP <- plot_half_violin(
  dat_MAP, "Blood", "IGFL3", "Clinical_diagnosis", unit_label = " (MAP)"
)
