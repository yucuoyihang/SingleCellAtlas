# Cell-fate signature scoring and visualization.
# Finalized gene sets are scored on the rooted Monocle2 CellDataSet with AUCell.

suppressPackageStartupMessages({
  library(monocle)
  library(AUCell)
  library(GSEABase)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(grid)
})


validate_signatures <- function(signatures) {
  signatures <- copy(as.data.table(signatures))
  required_cols <- c("signature_id", "signature_label", "gene_symbol")
  missing_cols <- setdiff(required_cols, colnames(signatures))
  if (length(missing_cols) > 0) stop("signature table missing: ", paste(missing_cols, collapse = ", "))
  for (column in required_cols) {
    set(signatures, j = column, value = as.character(signatures[[column]]))
    if (anyNA(signatures[[column]]) || any(!nzchar(trimws(signatures[[column]])))) {
      stop("Signature identifiers, labels, and gene symbols must be nonmissing and nonempty")
    }
  }
  if (nrow(signatures) == 0L) stop("The finalized signature table is empty")
  signatures[, gene_symbol := toupper(trimws(gene_symbol))]
  if (!"category" %in% colnames(signatures)) signatures[, category := signature_label]
  labels <- unique(signatures[, .(signature_id, signature_label, category)])
  if (anyDuplicated(labels$signature_id) || anyDuplicated(labels$signature_label)) {
    stop("Each signature must have one distinct label and category")
  }
  unique(signatures)
}


extract_uppercase_expression <- function(cds, symbol_column = NULL) {
  expr <- exprs(cds)
  symbols <- rownames(expr)
  features <- fData(cds)
  if (is.null(symbols) || anyNA(symbols) || any(!nzchar(trimws(symbols)))) {
    stop("Expression features require nonmissing, nonempty identifiers")
  }
  if (is.null(symbol_column) && any(grepl("^ENS[A-Z]*G[0-9]+(\\.[0-9]+)?$", symbols))) {
    symbol_column <- "gene_short_name"
  }
  if (!is.null(symbol_column)) {
    if (!symbol_column %in% colnames(features) ||
        anyDuplicated(rownames(features)) || anyNA(match(symbols, rownames(features)))) {
      stop("Feature identifiers require an aligned fData symbol column: ", symbol_column)
    }
    symbols <- as.character(features[match(symbols, rownames(features)), symbol_column])
  }
  symbols <- toupper(trimws(symbols))
  if (anyNA(symbols) || any(!nzchar(symbols)) || any(grepl("^ENS[A-Z]*G[0-9]+(\\.[0-9]+)?$", symbols))) {
    stop("Expression gene symbols must be nonmissing and nonempty; map feature IDs using symbol_column")
  }
  if (!inherits(expr, "dgCMatrix")) expr <- as(expr, "CsparseMatrix")
  if (any(!is.finite(expr@x)) || any(expr@x < 0)) stop("Expression counts must be finite and nonnegative")
  # Features sharing an uppercase gene symbol contribute their summed counts.
  symbol_levels <- unique(symbols)
  membership <- sparseMatrix(i = match(symbols, symbol_levels), j = seq_along(symbols),
                             x = 1, dims = c(length(symbol_levels), length(symbols)))
  expr <- membership %*% expr
  rownames(expr) <- symbol_levels
  expr
}


validate_states <- function(states) {
  states <- as.character(states)
  if (length(states) == 0L || anyNA(states) || any(!states %in% c("1", "2", "3"))) {
    stop("State must contain only nonmissing values 1, 2, or 3 (CellFate1, CellFate2, Root)")
  }
  states
}


prepare_gene_sets <- function(signatures, expression_genes, minimum_overlap) {
  # At least five measured genes are required to score a signature by default.
  overlap <- signatures[, .(
    n_signature_genes = uniqueN(gene_symbol),
    n_overlap_genes = uniqueN(gene_symbol[gene_symbol %in% expression_genes]),
    overlap_genes = paste(sort(unique(gene_symbol[gene_symbol %in% expression_genes])), collapse = ";")
  ), by = .(signature_id, signature_label, category)]
  overlap[, keep_for_auc := n_overlap_genes >= minimum_overlap]
  setorder(overlap, -n_overlap_genes, signature_id)

  keep_ids <- overlap[keep_for_auc == TRUE, signature_id]
  if (length(keep_ids) == 0) stop("no signature passed the overlap filter")

  gene_sets <- unique(signatures[
    signature_id %in% keep_ids & gene_symbol %in% expression_genes,
    .(signature_id, gene_symbol)
  ])
  gene_set_list <- split(gene_sets$gene_symbol, gene_sets$signature_id)
  gene_set_list <- lapply(gene_set_list, unique)

  gene_sets_gsc <- GeneSetCollection(lapply(names(gene_set_list), function(nm) {
    GeneSet(gene_set_list[[nm]], setName = nm)
  }))
  list(overlap = overlap, gene_sets = gene_sets_gsc)
}


score_signatures_with_aucell <- function(expr, gene_sets, auc_rank_fraction = 0.05, seed = 123L) {
  # A fixed seed controls random ordering of tied expression values.
  set.seed(seed)
  rankings <- AUCell_buildRankings(
    expr,
    BPPARAM = BiocParallel::SerialParam(RNGseed = seed),
    plotStats = FALSE,
    splitByBlocks = TRUE,
    verbose = TRUE
  )
  # AUCell integrates recovery within the top 5% of ranked genes by default.
  auc_max_rank <- ceiling(auc_rank_fraction * nrow(expr))
  auc <- AUCell_calcAUC(
    gene_sets,
    rankings,
    nCores = 1,
    aucMaxRank = auc_max_rank,
    normAUC = TRUE,
    verbose = TRUE
  )
  list(auc = auc, auc_max_rank = auc_max_rank)
}


merge_scores_with_metadata <- function(cds, auc) {
  auc_matrix <- as.data.frame(t(getAUC(auc)))
  auc_dt <- as.data.table(auc_matrix, keep.rownames = "cell")
  metadata <- as.data.table(pData(cds), keep.rownames = "cell")
  merge(metadata, auc_dt, by = "cell", all.x = TRUE, sort = FALSE)
}


summarize_signature_scores <- function(score_meta, score_cols, signatures) {
  id_vars <- intersect(c("cell", "fine", "origin", "b_diagnosis", "c_diagnosis", "State", "Pseudotime"), colnames(score_meta))
  score_long <- melt(
    score_meta[, c(id_vars, score_cols), with = FALSE],
    id.vars = id_vars,
    measure.vars = score_cols,
    variable.name = "signature_id",
    value.name = "auc_score"
  )
  labels <- signatures[, .(
    signature_label = signature_label[1],
    category = paste(sort(unique(category)), collapse = "; ")
  ), by = signature_id]
  score_long <- merge(score_long, labels, by = "signature_id", all.x = TRUE)

  summaries <- list()
  for (grp in intersect(c("fine", "origin", "b_diagnosis", "c_diagnosis", "State"), colnames(score_long))) {
    tmp <- score_long[!is.na(get(grp)), .(
      n_cells = .N,
      mean_auc = mean(auc_score, na.rm = TRUE),
      median_auc = median(auc_score, na.rm = TRUE),
      sd_auc = sd(auc_score, na.rm = TRUE)
    ), by = c(grp, "signature_id", "signature_label", "category")]
    setnames(tmp, grp, "group_value")
    tmp[, group_var := grp]
    summaries[[grp]] <- tmp
  }
  if ("Pseudotime" %in% colnames(score_long)) {
    score_long[, pseudotime_bin := cut(Pseudotime, breaks = 20, include.lowest = TRUE, ordered_result = TRUE)]
    pt <- score_long[!is.na(pseudotime_bin), .(
      n_cells = .N,
      mean_auc = mean(auc_score, na.rm = TRUE),
      median_auc = median(auc_score, na.rm = TRUE),
      sd_auc = sd(auc_score, na.rm = TRUE),
      pseudotime_mid = mean(Pseudotime, na.rm = TRUE)
    ), by = .(pseudotime_bin, signature_id, signature_label, category)]
    pt[, `:=`(group_var = "pseudotime_bin", group_value = as.character(pseudotime_bin))]
    summaries[["pseudotime_bin"]] <- pt[, .(group_var, group_value, signature_id, signature_label, category, n_cells, mean_auc, median_auc, sd_auc, pseudotime_mid)]
  }
  list(long = score_long, summary = rbindlist(summaries, fill = TRUE))
}


signature_function_classes <- function() {
  data.table(
    signature_id = c(
      "Cellular_senescence", "SASP", "T_cell_cytotoxicity", "Leukocyte_cytotoxicity",
      "Integrin_adhesion", "Cell_matrix_adhesion", "Leukocyte_endothelial_adhesion",
      "T_cell_migration", "TCR_signaling", "BBB_integrity", "Tight_junction_organization",
      "BBB_permeability_opening", "Neuroinflammatory_response", "Interferon_gamma_response",
      "Complement", "Phagocytosis", "Lysosome", "Amyloid_beta_clearance",
      "Oxidative_stress", "Oxidative_phosphorylation", "Glycolysis", "Apoptosis",
      "T_cell_exhaustion"
    ),
    function_class = c(
      "Senescence/Cell death", "Senescence/Cell death", "Immune signaling", "Immune signaling",
      "Barrier/Adhesion", "Barrier/Adhesion", "Immune signaling", "Immune signaling",
      "Immune signaling", "Barrier/Adhesion", "Barrier/Adhesion", "Barrier/Adhesion",
      "Immune signaling", "Immune signaling", "Immune signaling", "Clearance", "Clearance",
      "Clearance", "Metabolism/Stress", "Metabolism/Stress", "Metabolism/Stress",
      "Senescence/Cell death", "Immune signaling"
    )
  )
}


safe_group_zscore <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sigma <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sigma) || sigma == 0) return(rep(0, length(x)))
  (x - mu) / sigma
}


build_dotplot_tables <- function(score_meta, overlap) {
  score_meta <- copy(score_meta)
  state_label_map <- c("3" = "Root", "1" = "CellFate1", "2" = "CellFate2")
  score_cols <- overlap[keep_for_auc == TRUE][order(-n_overlap_genes), signature_id]
  score_cols <- score_cols[score_cols %in% colnames(score_meta)]

  label_map <- overlap[, .(
    signature_label = signature_label[1],
    n_overlap_genes = max(n_overlap_genes)
  ), by = signature_id]
  label_map <- merge(label_map, signature_function_classes(), by = "signature_id", all.x = TRUE, sort = FALSE)
  label_map[is.na(function_class), function_class := "Other"]

  score_meta[, State_label := state_label_map[validate_states(State)]]
  score_meta[, State_label := factor(State_label, levels = unname(state_label_map))]

  long <- melt(
    score_meta[, c("cell", "State_label", "fine", score_cols), with = FALSE],
    id.vars = c("cell", "State_label", "fine"),
    measure.vars = score_cols,
    variable.name = "signature_id",
    value.name = "auc_score"
  )
  long <- merge(long, label_map, by = "signature_id", all.x = TRUE)
  # Point size is the percentage above the pooled per-signature median,
  # not the percentage exceeding a biological activity threshold.
  long[, global_median := median(auc_score, na.rm = TRUE), by = signature_id]

  state_dt <- long[, .(
    mean_auc = mean(auc_score, na.rm = TRUE),
    pct_above_median = mean(auc_score > global_median, na.rm = TRUE) * 100,
    n_cells = .N
  ), by = .(State_label, signature_id, signature_label, function_class)]
  # Per-signature z scores compare state means, or state-by-fine means below;
  # these are relative scores, clipped to [-2.5, 2.5] for color display.
  state_dt[, z_mean_auc := safe_group_zscore(mean_auc), by = signature_id]
  state_dt[, z_mean_auc := pmax(pmin(z_mean_auc, 2.5), -2.5)]

  state_fine_dt <- long[, .(
    mean_auc = mean(auc_score, na.rm = TRUE),
    pct_above_median = mean(auc_score > global_median, na.rm = TRUE) * 100,
    n_cells = .N
  ), by = .(State_label, fine, signature_id, signature_label, function_class)]
  state_fine_dt[, z_mean_auc := safe_group_zscore(mean_auc), by = signature_id]
  state_fine_dt[, z_mean_auc := pmax(pmin(z_mean_auc, 2.5), -2.5)]

  cluster_dt <- dcast(state_dt, signature_id ~ State_label, value.var = "z_mean_auc", fill = 0)
  cluster_mat <- as.matrix(cluster_dt[, -1, with = FALSE])
  rownames(cluster_mat) <- cluster_dt$signature_id
  cluster_mat[!is.finite(cluster_mat)] <- 0
  signature_order_ids <- if (nrow(cluster_mat) > 1) {
    rownames(cluster_mat)[hclust(dist(cluster_mat), method = "average")$order]
  } else {
    rownames(cluster_mat)
  }
  signature_levels <- rev(label_map[match(signature_order_ids, signature_id), signature_label])
  list(state = state_dt, state_fine = state_fine_dt, signature_levels = signature_levels)
}


display_fine_label <- function(x) {
  as.character(x)
}


plot_signature_dotplot <- function(dot_tables, level = c("state", "state_fine")) {
  level <- match.arg(level)
  by_fine <- level == "state_fine"
  class_levels <- c("Immune signaling", "Barrier/Adhesion", "Clearance", "Metabolism/Stress", "Senescence/Cell death", "Other")
  colors <- if (by_fine) {
    c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")
  } else {
    c("#001A33", "#084594", "#FFFFFF", "#B2182B", "#67001F")
  }
  dt <- copy(dot_tables[[level]])
  dt[, State_label := factor(State_label, levels = c("Root", "CellFate1", "CellFate2"))]
  dt[, signature_label := factor(signature_label, levels = dot_tables$signature_levels)]
  dt[, function_class_plot := factor(toupper(function_class), levels = toupper(class_levels))]
  if (by_fine) dt[, fine := factor(fine, levels = sort(unique(as.character(fine))))]
  x_column <- if (by_fine) "fine" else "State_label"

  p <- ggplot(dt, aes(x = .data[[x_column]], y = signature_label)) +
    geom_point(aes(size = pct_above_median, fill = z_mean_auc), shape = 21,
               color = "#3F3F3F", stroke = if (by_fine) 0.24 else 0.28) +
    facet_grid(if (by_fine) function_class_plot ~ State_label else function_class_plot ~ .,
               scales = if (by_fine) "free" else "free_y",
               space = if (by_fine) "free" else "free_y") +
    scale_fill_gradientn(colors = colors, limits = c(-2.5, 2.5),
                         breaks = c(-2, -1, 0, 1, 2), name = "Mean score\nz-score") +
    scale_size_continuous(range = if (by_fine) c(0.8, 5.5) else c(1.2, 6.5),
                          limits = c(0, 100), name = "% cells above\nsignature median") +
    guides(
      fill = guide_colorbar(order = 1, title.position = "top", barheight = unit(3.8, "cm"), barwidth = unit(0.45, "cm")),
      size = guide_legend(order = 2, title.position = "top", override.aes = list(shape = 21, fill = "#FFFFFF", color = "black", stroke = 0.8))
    ) +
    theme_bw(base_family = "Arial", base_size = if (by_fine) 9 else 10) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(family = "Arial", color = "black", size = 10),
      strip.text.y = element_text(family = "Arial", color = "black", angle = 0, size = 10),
      panel.grid.major = element_line(color = "grey82", linetype = "dashed", linewidth = if (by_fine) 0.22 else 0.25),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = if (by_fine) 55 else 35, hjust = 1, vjust = 1,
                                 family = "Arial", color = "black", size = if (by_fine) 9 else 10),
      axis.text.y = element_text(family = "Arial", color = "black"),
      axis.title = element_blank(),
      plot.title = if (by_fine) element_text() else element_blank(),
      legend.title = element_text(family = "Arial", face = "bold", size = 11),
      legend.text = element_text(family = "Arial", size = 10),
      legend.position = "right",
      panel.border = element_rect(color = "black", fill = NA, linewidth = if (by_fine) 0.55 else 0.6),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  if (by_fine) {
    p <- p + scale_x_discrete(labels = display_fine_label) +
      labs(title = "Signature activity by fine annotation within root state and cell fates")
  }
  p
}


# Inputs are a rooted CellDataSet and a finalized table of signature_id,
# signature_label, and gene_symbol; category and source annotations are optional.
# Non-Ensembl feature identifiers can be mapped with an explicit symbol_column.
run_signature_analysis <- function(cds, signatures, minimum_overlap_genes = 5L,
                                   auc_rank_fraction = 0.05, seed = 123L,
                                   symbol_column = NULL, source_annotations = NULL) {
  required_metadata <- c("fine", "Pseudotime", "State")
  missing_metadata <- setdiff(required_metadata, colnames(pData(cds)))
  if (length(missing_metadata) > 0) stop("CDS missing metadata: ", paste(missing_metadata, collapse = ", "))
  validate_states(pData(cds)$State)
  if (length(minimum_overlap_genes) != 1L || !is.finite(minimum_overlap_genes) ||
      minimum_overlap_genes < 1 || minimum_overlap_genes != floor(minimum_overlap_genes)) {
    stop("minimum_overlap_genes must be a positive integer")
  }
  if (length(auc_rank_fraction) != 1L || !is.finite(auc_rank_fraction) ||
      auc_rank_fraction <= 0 || auc_rank_fraction > 1) {
    stop("auc_rank_fraction must be in (0, 1]")
  }
  if (length(seed) != 1L || !is.finite(seed) || seed < 0 ||
      seed > .Machine$integer.max || seed != floor(seed)) stop("seed must be a nonnegative integer")
  signatures <- validate_signatures(signatures)
  if (any(signatures$signature_id %in% c(colnames(pData(cds)), "cell", "State_label"))) {
    stop("Signature identifiers must not collide with cell metadata columns")
  }
  expr <- extract_uppercase_expression(cds, symbol_column)
  prepared <- prepare_gene_sets(signatures, rownames(expr), minimum_overlap_genes)
  auc_result <- score_signatures_with_aucell(expr, prepared$gene_sets, auc_rank_fraction, seed)
  score_meta <- merge_scores_with_metadata(cds, auc_result$auc)
  summaries <- summarize_signature_scores(score_meta, rownames(getAUC(auc_result$auc)), signatures)
  dot_tables <- build_dotplot_tables(score_meta, prepared$overlap)

  result <- list(
    signatures = signatures,
    gene_sets = prepared$gene_sets,
    overlap = prepared$overlap,
    auc = auc_result$auc,
    auc_max_rank = auc_result$auc_max_rank,
    score_metadata = score_meta,
    summaries = summaries,
    dot_tables = dot_tables,
    plots = list(state = plot_signature_dotplot(dot_tables, "state"),
                 state_fine = plot_signature_dotplot(dot_tables, "state_fine"))
  )
  if (!is.null(source_annotations)) attr(result, "source_annotations") <- source_annotations
  result
}
