# Monocle2 trajectory inference and branch-dependent expression analysis.
# Ordering genes, DDRTree, and BEAM follow the Monocle2 analysis framework.
# State identifiers describe inferred graph segments, not validated cell fates.

library(monocle)
library(Biobase)
library(Matrix)
library(data.table)
library(ggplot2)

validate_monocle2_inputs <- function(expr_matrix, metadata, gene_metadata) {
  if (!(is.matrix(expr_matrix) || inherits(expr_matrix, "Matrix")) ||
      length(dim(expr_matrix)) != 2L || any(dim(expr_matrix) == 0L)) {
    stop("expr_matrix must be a nonempty gene-by-cell count matrix.")
  }
  metadata <- as.data.frame(metadata)
  gene_metadata <- as.data.frame(gene_metadata)
  if (!all(c("barcode", "fine") %in% names(metadata))) {
    stop("metadata must contain barcode and fine.")
  }
  if (!"gene_id" %in% names(gene_metadata)) {
    stop("gene_metadata must contain gene_id.")
  }
  valid_ids <- function(x) {
    !anyNA(x) && all(nzchar(as.character(x))) && !anyDuplicated(x)
  }
  if (!valid_ids(metadata$barcode) || !valid_ids(gene_metadata$gene_id)) {
    stop("Metadata barcode and gene_id values must be nonmissing, nonempty, and unique.")
  }
  metadata$barcode <- as.character(metadata$barcode)
  gene_metadata$gene_id <- as.character(gene_metadata$gene_id)
  if (anyNA(metadata$fine) || any(!nzchar(as.character(metadata$fine)))) {
    stop("fine must contain nonmissing, nonempty sampling strata.")
  }
  if (nrow(expr_matrix) != nrow(gene_metadata) ||
      ncol(expr_matrix) != nrow(metadata)) {
    stop("Matrix dimensions must match the gene and cell metadata.")
  }
  # Named matrix axes determine metadata alignment; unnamed axes use metadata order.
  if (!is.null(rownames(expr_matrix))) {
    if (!valid_ids(rownames(expr_matrix)) ||
        !setequal(rownames(expr_matrix), gene_metadata$gene_id)) {
      stop("Matrix row names must uniquely match gene_metadata$gene_id.")
    }
    gene_metadata <- gene_metadata[match(rownames(expr_matrix), gene_metadata$gene_id), , drop = FALSE]
  }
  if (!is.null(colnames(expr_matrix))) {
    if (!valid_ids(colnames(expr_matrix)) ||
        !setequal(colnames(expr_matrix), metadata$barcode)) {
      stop("Matrix column names must uniquely match metadata$barcode.")
    }
    metadata <- metadata[match(colnames(expr_matrix), metadata$barcode), , drop = FALSE]
  }
  if (!"gene_short_name" %in% names(gene_metadata)) {
    gene_metadata$gene_short_name <- gene_metadata$gene_id
  }
  gene_metadata$gene_short_name <- as.character(gene_metadata$gene_short_name)
  missing_names <- is.na(gene_metadata$gene_short_name) | !nzchar(gene_metadata$gene_short_name)
  gene_metadata$gene_short_name[missing_names] <- gene_metadata$gene_id[missing_names]
  rownames(expr_matrix) <- rownames(gene_metadata) <- gene_metadata$gene_id
  colnames(expr_matrix) <- rownames(metadata) <- metadata$barcode
  values <- if (inherits(expr_matrix, "sparseMatrix")) {
    if (!"x" %in% slotNames(expr_matrix)) stop("expr_matrix must contain numeric counts.")
    expr_matrix@x
  } else {
    as.vector(expr_matrix)
  }
  if (!is.numeric(values) || any(!is.finite(values)) || any(values < 0) ||
      any(abs(values - round(values)) > 1e-8)) {
    stop("expr_matrix must contain finite, nonnegative integer counts.")
  }
  if (any(Matrix::colSums(expr_matrix) == 0)) stop("Each cell must have a positive count total.")
  list(expr_matrix = expr_matrix, metadata = metadata, gene_metadata = gene_metadata)
}

validate_scalar <- function(x, name, lower, upper = Inf, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < lower || x > upper || (integer && x != floor(x))) {
    stop(name, " must be a finite ", if (integer) "integer " else "",
         "value in [", lower, ", ", upper, "].")
  }
}

stratified_sample_by_fine <- function(expr_matrix, metadata, fraction, seed) {
  validate_scalar(fraction, "sample_fraction", .Machine$double.eps, 1)
  validate_scalar(seed, "random_seed", 0, .Machine$integer.max, integer = TRUE)
  # Each fine annotation retains 30% of its cells by default, with at least one cell per stratum.
  if (fraction == 1) return(list(expr_matrix = expr_matrix, metadata = metadata))
  set.seed(seed)
  strata <- split(seq_len(nrow(metadata)),
                  factor(metadata$fine, levels = unique(metadata$fine)))
  selected <- unlist(lapply(strata, function(i) {
    i[sample.int(length(i), max(1L, round(length(i) * fraction)))]
  }), use.names = FALSE)
  selected <- sort(selected)
  list(expr_matrix = expr_matrix[, selected, drop = FALSE],
       metadata = metadata[selected, , drop = FALSE])
}

build_monocle2_cds <- function(expr_matrix, metadata, gene_metadata) {
  cds <- newCellDataSet(
    as(expr_matrix, "sparseMatrix"),
    phenoData = new("AnnotatedDataFrame", data = as.data.frame(metadata)),
    featureData = new("AnnotatedDataFrame", data = as.data.frame(gene_metadata)),
    expressionFamily = negbinomial.size()
  )
  detectGenes(cds, min_expr = 0.1)
}

run_monocle2_ordering <- function(cds, random_seed) {
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  ordering_genes <- rownames(fData(cds))
  if (length(ordering_genes) != 3000L) stop("The CDS must contain exactly 3,000 highly variable genes.")
  cds <- setOrderingFilter(cds, ordering_genes)
  set.seed(random_seed)
  cds <- reduceDimension(cds, max_components = 2, method = "DDRTree", verbose = TRUE)
  orderCells(cds)
}

validate_branch_selection <- function(states, root_state, branch_point, branch_states,
                                      n_branch_points) {
  if (length(root_state) != 1L || anyNA(root_state) ||
      !as.character(root_state) %in% as.character(states)) {
    stop("root_state must be present in the ordered CDS; available states: ",
         paste(states, collapse = ", "), ".")
  }
  validate_scalar(branch_point, "branch_point", 1, integer = TRUE)
  if (branch_point > n_branch_points) {
    stop("branch_point ", branch_point, " is absent; the DDRTree graph has ",
         n_branch_points, " branch points.")
  }
  if (!is.null(branch_states)) {
    if (length(branch_states) != 2L || anyNA(branch_states) ||
        anyDuplicated(as.character(branch_states)) ||
        !all(as.character(branch_states) %in% as.character(states))) {
      stop("branch_states must contain two distinct available states; available states: ",
           paste(states, collapse = ", "), ".")
    }
    if (as.character(root_state) %in% as.character(branch_states)) {
      stop("branch_states must exclude the root state.")
    }
  }
}

run_rooted_beam <- function(cds, root_state, branch_point, branch_states, branch_labels) {
  validate_branch_selection(unique(pData(cds)$State), root_state, branch_point,
                            branch_states = NULL,
                            n_branch_points = length(cds@auxOrderingData[["DDRTree"]]$branch_points))
  cds_root <- orderCells(cds, root_state = root_state)
  rooted_pd <- pData(cds_root)
  if (any(!is.finite(rooted_pd$Pseudotime))) stop("Rooted pseudotime contains nonfinite values.")
  actual_root <- unique(as.character(rooted_pd$State[rooted_pd$Pseudotime == 0]))
  if (length(actual_root) != 1L) stop("Rooted CDS must identify a single state at pseudotime zero.")
  validate_branch_selection(unique(rooted_pd$State), actual_root, branch_point, branch_states,
                            length(cds_root@auxOrderingData[["DDRTree"]]$branch_points))
  # Explicit state selection takes precedence over branch_point in Monocle 2.26.0.
  # Duplicated progenitor cells match the branched-heatmap model specification.
  branch_cds <- tryCatch(
    buildBranchCellDataSet(cds_root, branch_point = branch_point,
                          branch_states = branch_states, branch_labels = branch_labels,
                          progenitor_method = "duplicate"),
    error = function(e) stop("Cannot construct the selected root-to-branch paths: ",
                             conditionMessage(e), call. = FALSE)
  )
  branch_pd <- pData(branch_cds)
  if (length(unique(branch_pd$Branch)) != 2L ||
      any(!is.finite(branch_pd$Pseudotime))) {
    stop("Selected paths must yield two branches with finite pseudotime.")
  }
  beam_res <- BEAM(cds_root, branch_point = branch_point, branch_states = branch_states,
                   branch_labels = branch_labels, progenitor_method = "duplicate",
                   cores = 1, verbose = TRUE)
  beam_dt <- as.data.table(beam_res)
  if (!"gene_id" %in% names(beam_dt)) beam_dt[, gene_id := rownames(beam_res)]
  if (!"gene_short_name" %in% names(beam_dt)) {
    beam_dt[, gene_short_name := fData(cds_root)[gene_id, "gene_short_name"]]
  }
  setorder(beam_dt, qval, pval)
  list(cds_root = cds_root, beam = beam_dt)
}

select_beam_heatmap_genes <- function(beam_dt, qval_threshold, max_genes) {
  validate_scalar(qval_threshold, "beam_qval_threshold", 0, 1)
  validate_scalar(max_genes, "max_heatmap_genes", 2, integer = TRUE)
  selected <- beam_dt[is.finite(qval) & qval < qval_threshold][order(qval, pval)]
  if (nrow(selected) == 0L) stop("No BEAM genes passed the specified q-value threshold.")
  selected$gene_id[seq_len(min(nrow(selected), max_genes))]
}

plot_monocle2_trajectory <- function(cds, color_by = "fine") {
  if (!color_by %in% colnames(pData(cds))) stop("Trajectory color_by is absent from cell metadata: ", color_by)
  plot_cell_trajectory(cds, color_by = color_by, show_branch_points = TRUE) +
    theme_classic(base_family = "Arial") +
    theme(axis.title = element_text(family = "Arial"),
          axis.text = element_text(family = "Arial", color = "black"),
          legend.title = element_text(family = "Arial"),
          legend.text = element_text(family = "Arial"),
          plot.title = element_blank())
}

calculate_beam_heatmap <- function(cds_root, heatmap_genes, branch_point,
                                   branch_states, branch_labels, heatmap_clusters) {
  if (length(heatmap_genes) < max(2L, heatmap_clusters)) {
    stop("Too few selected BEAM genes for the requested heatmap clustering.")
  }
  # Stable gene identifiers preserve distinct rows when gene symbols are duplicated.
  # The native routine duplicates progenitors and draws its return object; a null
  # graphics device retains the calculation without displaying an intermediate figure.
  grDevices::pdf(file = NULL)
  device <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(device), add = TRUE)
  tryCatch(
    plot_genes_branched_heatmap(
      cds_root[heatmap_genes, ], branch_point = branch_point,
      branch_states = branch_states, branch_labels = branch_labels,
      num_clusters = heatmap_clusters, cores = 1, show_rownames = FALSE,
      use_gene_short_name = FALSE, return_heatmap = TRUE,
      branch_colors = c("#797D7D", "#A71D1D", "#234862"),
      hmcols = grDevices::colorRampPalette(c("#234862", "white", "#A71D1D"))(100)
    ),
    error = function(e) stop("Branched heatmap calculation failed; check fitted expression and ",
                             "the number of nonconstant genes for clustering: ",
                             conditionMessage(e), call. = FALSE)
  )
}

extract_heatmap_clusters <- function(hm, beam_dt) {
  cluster_dt <- data.table(
    gene_id = rownames(hm$annotation_row),
    Cluster = paste0("Cluster_", as.character(hm$annotation_row$Cluster))
  )
  merge(cluster_dt, beam_dt, by = "gene_id", all.x = TRUE, sort = FALSE)
}

run_monocle2_enrichment <- function(heatmap_clusters, beam_dt,
                                    pvalue_cutoff = 0.05,
                                    min_gene_set_size = 10L,
                                    max_gene_set_size = 500L) {
  # Module assignments originate from the BEAM branched-heatmap clustering.
  # gene_short_name contains human gene symbols, independently of gene_id format.
  heatmap_clusters <- as.data.frame(heatmap_clusters)
  beam_dt <- as.data.frame(beam_dt)
  stopifnot(all(c("gene_id", "Cluster") %in% names(heatmap_clusters)),
            all(c("gene_id", "gene_short_name", "pval") %in% names(beam_dt)))
  validate_scalar(pvalue_cutoff, "pvalue_cutoff", 0, 1)
  validate_scalar(min_gene_set_size, "min_gene_set_size", 1, integer = TRUE)
  validate_scalar(max_gene_set_size, "max_gene_set_size", min_gene_set_size, integer = TRUE)
  if (!nrow(heatmap_clusters) || anyNA(heatmap_clusters$Cluster) ||
      any(!nzchar(as.character(heatmap_clusters$Cluster)))) {
    stop("BEAM modules must have nonmissing, nonempty cluster labels.")
  }

  # The universe comprises successfully tested genes, not only significant BEAM genes.
  tested <- is.finite(beam_dt$pval)
  if ("status" %in% names(beam_dt)) {
    tested <- tested & !is.na(beam_dt$status) & beam_dt$status == "OK"
  }
  gene_metadata <- unique(beam_dt[tested, c("gene_id", "gene_short_name")])
  gene_metadata <- gene_metadata[!is.na(gene_metadata$gene_short_name) &
                                   nzchar(as.character(gene_metadata$gene_short_name)), , drop = FALSE]
  if (!nrow(gene_metadata)) stop("No tested genes with human gene symbols are available.")
  gene_mapping <- clusterProfiler::bitr(
    unique(as.character(gene_metadata$gene_short_name)),
    fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )
  gene_mapping <- merge(gene_metadata, gene_mapping,
                         by.x = "gene_short_name", by.y = "SYMBOL")
  universe <- unique(as.character(gene_mapping$ENTREZID))
  if (!length(universe)) stop("No tested genes mapped to Entrez identifiers.")

  # Retain all valid mappings and count each Entrez identifier once per module.
  module_genes <- split(as.character(heatmap_clusters$gene_id),
                         as.character(heatmap_clusters$Cluster))
  module_genes <- lapply(module_genes, function(gene_ids) {
    unique(as.character(gene_mapping$ENTREZID[gene_mapping$gene_id %in% gene_ids]))
  })
  # BH correction is applied separately within each module and annotation database.
  # qvalueCutoff = 1 avoids an additional q-value threshold beyond the specified P cutoff.
  go_bp <- lapply(module_genes, function(genes) {
    if (!length(genes)) return(NULL)
    clusterProfiler::enrichGO(
      gene = genes, universe = universe, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
      pvalueCutoff = pvalue_cutoff, qvalueCutoff = 1,
      minGSSize = min_gene_set_size, maxGSSize = max_gene_set_size, readable = TRUE
    )
  })
  # KEGG annotations are retrieved online for Homo sapiens.
  kegg <- lapply(module_genes, function(genes) {
    if (!length(genes)) return(NULL)
    clusterProfiler::enrichKEGG(
      gene = genes, universe = universe, organism = "hsa", keyType = "ncbi-geneid",
      pAdjustMethod = "BH", pvalueCutoff = pvalue_cutoff, qvalueCutoff = 1,
      minGSSize = min_gene_set_size, maxGSSize = max_gene_set_size,
      use_internal_data = FALSE
    )
  })
  list(go_bp = go_bp, kegg = kegg, module_genes = module_genes,
       universe = universe, gene_mapping = gene_mapping)
}

plot_protein_labeled_heatmap <- function(hm, gene_metadata, protein_genes, heatmap_clusters) {
  # Protein-associated symbols label matching rows without removing other BEAM genes.
  symbols <- gene_metadata$gene_short_name[match(rownames(hm$heatmap_matrix), gene_metadata$gene_id)]
  label_genes <- intersect(protein_genes, symbols)
  if (length(label_genes) == 0L) stop("No protein_genes symbols matched the retained BEAM heatmap rows.")
  row_labels <- ifelse(symbols %in% label_genes, symbols, "")
  # The native heatmap uses 0.1-spaced breaks; the labeled view retains this scale.
  expression_range <- range(hm$heatmap_matrix)
  pheatmap::pheatmap(
    hm$heatmap_matrix, color = hm$hmcols,
    breaks = seq(expression_range[1] - 0.1, expression_range[2] + 0.1, by = 0.1),
    cluster_rows = hm$ph_res$tree_row, cutree_rows = heatmap_clusters,
    cluster_cols = FALSE, show_rownames = TRUE, labels_row = row_labels,
    show_colnames = FALSE, annotation_row = hm$annotation_row,
    annotation_col = hm$annotation_col, annotation_colors = hm$annotation_colors,
    gaps_col = hm$col_gap_ind, border_color = NA, fontsize_row = 6.4,
    main = "BEAM branched heatmap: protein-intersect genes labeled", silent = TRUE
  )
}

run_monocle2_analysis <- function(expr_matrix, metadata, gene_metadata, protein_genes,
                                  sample_fraction = 0.30, random_seed = 20260603,
                                  root_state = 3L, branch_point = 1L, branch_states = c(1, 2),
                                  beam_qval_threshold = 1e-4, max_heatmap_genes = 300L,
                                  heatmap_clusters = 5L) {
  # Default parameters retain State 3 rooting and States 1 and 2 as graph-defined branches.
  # BEAM genes satisfy q < 1e-4; at most 300 genes form five expression-pattern clusters.
  # protein_genes is a character vector of gene symbols from the supplied protein intersection.
  inputs <- validate_monocle2_inputs(expr_matrix, metadata, gene_metadata)
  validate_scalar(sample_fraction, "sample_fraction", .Machine$double.eps, 1)
  validate_scalar(random_seed, "random_seed", 0, .Machine$integer.max, integer = TRUE)
  validate_scalar(branch_point, "branch_point", 1, integer = TRUE)
  validate_scalar(beam_qval_threshold, "beam_qval_threshold", 0, 1)
  validate_scalar(max_heatmap_genes, "max_heatmap_genes", 2, integer = TRUE)
  validate_scalar(heatmap_clusters, "heatmap_clusters", 2, max_heatmap_genes, integer = TRUE)
  if (!is.character(protein_genes) || length(protein_genes) == 0L ||
      anyNA(protein_genes) || any(!nzchar(protein_genes))) {
    stop("protein_genes must be a nonempty character vector of nonmissing gene symbols.")
  }
  branch_labels <- if (is.null(branch_states)) c("Cell fate 1", "Cell fate 2") else paste0("State", branch_states, " fate")
  sampled <- stratified_sample_by_fine(inputs$expr_matrix, inputs$metadata,
                                       sample_fraction, random_seed)
  # VST ranks count variability on the sampled cells before restricting the CDS.
  # Temporary identifiers prevent Seurat feature-name sanitization from altering gene alignment.
  hvg_counts <- sampled$expr_matrix
  feature_ids <- paste0("feature", seq_len(nrow(hvg_counts)))
  rownames(hvg_counts) <- feature_ids
  assay <- SeuratObject::CreateAssayObject(counts = hvg_counts)
  assay <- Seurat::FindVariableFeatures(assay, selection.method = "vst",
                                        nfeatures = 3000L, verbose = FALSE)
  ordering_genes <- rownames(sampled$expr_matrix)[
    match(SeuratObject::VariableFeatures(assay), feature_ids)
  ]
  if (length(ordering_genes) != 3000L) {
    stop("The sampled expression matrix must yield exactly 3,000 highly variable genes.")
  }
  cds <- build_monocle2_cds(sampled$expr_matrix[ordering_genes, , drop = FALSE],
                            sampled$metadata,
                            inputs$gene_metadata[ordering_genes, , drop = FALSE])
  cds <- run_monocle2_ordering(cds, random_seed)
  beam_run <- run_rooted_beam(cds, root_state, branch_point, branch_states, branch_labels)
  heatmap_genes <- select_beam_heatmap_genes(beam_run$beam, beam_qval_threshold, max_heatmap_genes)
  hm <- calculate_beam_heatmap(beam_run$cds_root, heatmap_genes, branch_point,
                               branch_states, branch_labels, heatmap_clusters)
  list(
    cds_root = beam_run$cds_root,
    ordering_genes = rownames(fData(cds))[fData(cds)$use_for_ordering],
    sampled_cells = colnames(sampled$expr_matrix),
    beam = beam_run$beam, heatmap_genes = heatmap_genes,
    heatmap_clusters = extract_heatmap_clusters(hm, beam_run$beam),
    branched_heatmap = hm,
    trajectory_plot = plot_monocle2_trajectory(beam_run$cds_root, color_by = "fine"),
    protein_heatmap = plot_protein_labeled_heatmap(hm, inputs$gene_metadata,
                                                  unique(protein_genes), heatmap_clusters)
  )
}
