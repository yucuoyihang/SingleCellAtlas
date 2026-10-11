# DEG trajectory spectral clustering (CLS)

CLS_KEY_CELLS <- c(
  "CD8T_Teff_01_GZMH",
  "CD8T_Teff_02_KLRG1",
  "ncMono_04_BCL2",
  "ncMono_05_LYVE1"
)

AD_CONTINUUM_STAGES <- c(
  "CN A-",
  "Preclinical AD",
  "MCI due to AD",
  "Dementia due to AD"
)

ADJACENT_STAGE_COMPARISONS <- c(
  "Stage1_vs_Stage0",
  "Stage2_vs_Stage1",
  "Stage3_vs_Stage2"
)

require_cls_columns <- function(data, columns, object_name) {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      object_name,
      " is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

require_cls_packages <- function() {
  packages <- c("dplyr", "tidyr", "SNFtool")
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages) > 0) {
    stop(
      "Missing required packages: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }
}

# Select the union of significant genes from the three adjacent AD-continuum
# contrasts. pct.1 and pct.2 are expected to be proportions in [0, 1].
# Match Seurat min.pct semantics: expression in either group can satisfy the
# threshold, so use pmax() rather than requiring both groups with pmin().
select_adjacent_degs <- function(
    deg_results,
    log2fc_threshold = 0.25,
    adjusted_p_threshold = 0.05,
    min_pct = 0.10,
    adjacent_comparisons = ADJACENT_STAGE_COMPARISONS) {
  require_cls_columns(
    deg_results,
    c(
      "gene", "pair_comparison", "avg_log2FC", "p_val_adj",
      "pct.1", "pct.2"
    ),
    "deg_results"
  )

  numeric_columns <- c("avg_log2FC", "p_val_adj", "pct.1", "pct.2")
  non_numeric <- numeric_columns[
    !vapply(deg_results[numeric_columns], is.numeric, logical(1))
  ]
  if (length(non_numeric) > 0) {
    stop(
      "These DEG columns must be numeric: ",
      paste(non_numeric, collapse = ", "),
      call. = FALSE
    )
  }

  selected <- dplyr::filter(
    deg_results,
    .data$pair_comparison %in% adjacent_comparisons,
    abs(.data$avg_log2FC) > log2fc_threshold,
    .data$p_val_adj < adjusted_p_threshold,
    pmax(.data$pct.1, .data$pct.2) >= min_pct,
    !is.na(.data$gene)
  )

  genes <- sort(unique(as.character(selected$gene)))
  if (length(genes) == 0) {
    stop("No genes passed the adjacent-comparison DEG thresholds.", call. = FALSE)
  }

  selected <- dplyr::arrange(selected, .data$pair_comparison, .data$gene)
  list(rows = selected, genes = genes)
}

# Reproduce the trajectory feature construction used by the source analysis:
# Z-score each gene across samples, then average the Z-scores within each
# biological stage.
build_stage_profiles <- function(
    expression_long,
    sample_metadata,
    genes,
    stages = AD_CONTINUUM_STAGES) {
  require_cls_columns(
    expression_long,
    c("gene", "sample_id", "expression"),
    "expression_long"
  )
  require_cls_columns(
    sample_metadata,
    c("sample_id", "b_diagnosis"),
    "sample_metadata"
  )

  if (!is.numeric(expression_long$expression)) {
    stop("expression_long$expression must be numeric.", call. = FALSE)
  }
  if (anyDuplicated(expression_long[c("gene", "sample_id")])) {
    stop(
      "expression_long must contain one value per gene and sample.",
      call. = FALSE
    )
  }
  if (anyDuplicated(sample_metadata$sample_id)) {
    stop("sample_metadata must contain one row per sample_id.", call. = FALSE)
  }

  expression_input <- expression_long[
    expression_long$gene %in% genes &
      is.finite(expression_long$expression),
    c("gene", "sample_id", "expression"),
    drop = FALSE
  ]
  metadata_input <- sample_metadata[
    sample_metadata$b_diagnosis %in% stages,
    c("sample_id", "b_diagnosis"),
    drop = FALSE
  ]

  joined <- dplyr::inner_join(
    expression_input,
    metadata_input,
    by = "sample_id"
  )
  if (nrow(joined) == 0) {
    stop("No expression rows matched the requested genes and stages.", call. = FALSE)
  }

  available_stages <- unique(as.character(joined$b_diagnosis))
  missing_stages <- setdiff(stages, available_stages)
  if (length(missing_stages) > 0) {
    stop(
      "No matched samples for stages: ",
      paste(missing_stages, collapse = ", "),
      call. = FALSE
    )
  }

  stage_means <- joined |>
    dplyr::group_by(.data$gene) |>
    dplyr::mutate(
      gene_mean = mean(.data$expression),
      gene_sd = stats::sd(.data$expression),
      zscore = dplyr::if_else(
        is.finite(.data$gene_sd) & .data$gene_sd > 0,
        (.data$expression - .data$gene_mean) / .data$gene_sd,
        NA_real_
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(is.finite(.data$zscore)) |>
    dplyr::group_by(.data$gene, .data$b_diagnosis) |>
    dplyr::summarise(
      mean_zscore = mean(.data$zscore),
      .groups = "drop"
    )

  profiles <- tidyr::pivot_wider(
    stage_means,
    id_cols = "gene",
    names_from = "b_diagnosis",
    values_from = "mean_zscore"
  )
  profiles <- as.data.frame(profiles, check.names = FALSE)
  profiles <- profiles[, c("gene", stages), drop = FALSE]
  profiles <- profiles[stats::complete.cases(profiles), , drop = FALSE]

  if (nrow(profiles) == 0) {
    stop(
      "No genes had finite expression profiles across every stage.",
      call. = FALSE
    )
  }

  sample_counts <- joined |>
    dplyr::distinct(.data$sample_id, .data$b_diagnosis) |>
    dplyr::count(.data$b_diagnosis, name = "n_samples")
  sample_counts$b_diagnosis <- factor(
    sample_counts$b_diagnosis,
    levels = stages
  )
  sample_counts <- dplyr::arrange(sample_counts, .data$b_diagnosis)
  sample_counts$b_diagnosis <- as.character(sample_counts$b_diagnosis)

  list(
    profiles = profiles,
    sample_counts = sample_counts,
    n_selected_genes = length(unique(genes)),
    n_complete_genes = nrow(profiles)
  )
}

# Cluster genes by their stage-level Z-score trajectories. Euclidean
# dissimilarity is converted to the [0, 1] similarity used by SNFtool.
spectral_cluster_profiles <- function(
    stage_profiles,
    n_clusters = 5,
    stages = AD_CONTINUUM_STAGES,
    seed = 42) {
  require_cls_packages()
  require_cls_columns(
    stage_profiles,
    c("gene", stages),
    "stage_profiles"
  )

  if (length(n_clusters) != 1 || !is.finite(n_clusters) ||
      n_clusters < 2 || n_clusters != as.integer(n_clusters)) {
    stop("n_clusters must be one integer greater than 1.", call. = FALSE)
  }
  n_clusters <- as.integer(n_clusters)

  if (anyDuplicated(stage_profiles$gene)) {
    stop("stage_profiles$gene must be unique.", call. = FALSE)
  }
  if (nrow(stage_profiles) < n_clusters) {
    stop("The number of genes must be at least n_clusters.", call. = FALSE)
  }

  feature_matrix <- as.matrix(stage_profiles[, stages, drop = FALSE])
  storage.mode(feature_matrix) <- "double"
  rownames(feature_matrix) <- as.character(stage_profiles$gene)
  if (any(!is.finite(feature_matrix))) {
    stop("stage_profiles contains non-finite stage values.", call. = FALSE)
  }

  dissimilarity <- stats::dist(feature_matrix, method = "euclidean")
  max_distance <- max(as.numeric(dissimilarity))
  if (!is.finite(max_distance) || max_distance <= 0) {
    stop("All gene trajectories have zero pairwise distance.", call. = FALSE)
  }

  similarity <- 1 - as.matrix(dissimilarity) / max_distance
  similarity[similarity < 0] <- 0
  similarity[similarity > 1] <- 1
  diag(similarity) <- 1

  set.seed(seed)
  clusters <- SNFtool::spectralClustering(
    similarity,
    K = n_clusters
  )
  clusters <- as.integer(clusters)
  if (length(clusters) != nrow(feature_matrix) ||
      anyNA(clusters) ||
      any(!clusters %in% seq_len(n_clusters))) {
    stop("Spectral clustering returned invalid cluster labels.", call. = FALSE)
  }

  gene_profiles <- data.frame(
    gene = rownames(feature_matrix),
    cluster = clusters,
    feature_matrix,
    row.names = NULL,
    check.names = FALSE
  )
  gene_profiles <- dplyr::arrange(
    gene_profiles,
    .data$cluster,
    .data$gene
  )

  cluster_trajectories <- gene_profiles |>
    dplyr::group_by(.data$cluster) |>
    dplyr::summarise(
      n_genes = dplyr::n(),
      dplyr::across(dplyr::all_of(stages), mean),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$cluster)

  list(
    cluster_assignments = gene_profiles[, c("gene", "cluster")],
    gene_profiles = gene_profiles,
    cluster_trajectories = cluster_trajectories,
    similarity_matrix = similarity
  )
}

# End-to-end core analysis for one already prepared cell-type/tissue subset.
run_deg_cls_clustering <- function(
    deg_results,
    expression_long,
    sample_metadata,
    n_clusters = 5,
    log2fc_threshold = 0.25,
    adjusted_p_threshold = 0.05,
    min_pct = 0.10,
    stages = AD_CONTINUUM_STAGES,
    adjacent_comparisons = ADJACENT_STAGE_COMPARISONS,
    seed = 42) {
  require_cls_packages()

  selected <- select_adjacent_degs(
    deg_results = deg_results,
    log2fc_threshold = log2fc_threshold,
    adjusted_p_threshold = adjusted_p_threshold,
    min_pct = min_pct,
    adjacent_comparisons = adjacent_comparisons
  )
  stage_result <- build_stage_profiles(
    expression_long = expression_long,
    sample_metadata = sample_metadata,
    genes = selected$genes,
    stages = stages
  )
  cluster_result <- spectral_cluster_profiles(
    stage_profiles = stage_result$profiles,
    n_clusters = n_clusters,
    stages = stages,
    seed = seed
  )

  list(
    selected_degs = selected$rows,
    selected_genes = selected$genes,
    sample_counts = stage_result$sample_counts,
    stage_profiles = stage_result$profiles,
    cluster_assignments = cluster_result$cluster_assignments,
    gene_profiles = cluster_result$gene_profiles,
    cluster_trajectories = cluster_result$cluster_trajectories,
    similarity_matrix = cluster_result$similarity_matrix
  )
}
