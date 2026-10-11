# CellChat inference and SPP1 communication across tissue and diagnostic groups.
# Expression matrices are genes by cells; metadata row names are cell barcodes.

suppressPackageStartupMessages({
  library(Matrix)
  library(CellChat)
  library(data.table)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

signaling_of_interest <- "SPP1"
lyve1_sender <- "ncMono_05_LYVE1"


display_fine_label <- function(x) {
  as.character(x)
}


darken_hex <- function(cols, factor = 0.62) {
  vapply(cols, function(col) {
    if (is.na(col) || col == "") return("#555555")
    rgb <- grDevices::col2rgb(col)
    rgb <- pmax(0, pmin(255, round(rgb * factor)))
    grDevices::rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)
  }, character(1))
}


align_cellchat_metadata <- function(expr, metadata, columns = "fine") {
  metadata <- as.data.frame(metadata)
  if (!all(columns %in% names(metadata))) stop("Missing metadata columns: ", paste(setdiff(columns, names(metadata)), collapse = ", "))
  if (is.null(rownames(expr)) || anyNA(rownames(expr)) || anyDuplicated(rownames(expr)) || any(!nzchar(rownames(expr)))) {
    stop("Expression rows must have unique nonempty gene symbols.")
  }
  cells <- colnames(expr)
  if (is.null(cells) || anyNA(cells) || anyDuplicated(cells) || any(!nzchar(cells)) ||
      anyDuplicated(rownames(metadata)) || !setequal(cells, rownames(metadata))) {
    stop("Expression columns and metadata row names must contain the same unique barcodes.")
  }
  metadata <- metadata[cells, , drop = FALSE]
  if (any(vapply(metadata[columns], function(x) anyNA(x) || any(!nzchar(as.character(x))), logical(1)))) {
    stop("Grouping metadata cannot contain missing or empty values.")
  }
  metadata
}


run_cellchat_pipeline <- function(expr, metadata, database, group_col = "fine", min_cells = 10L) {
  metadata <- align_cellchat_metadata(expr, metadata, group_col)
  values <- if (inherits(expr, "sparseMatrix")) expr@x else as.vector(expr)
  if (any(!is.finite(values)) || any(values < 0)) stop("Expression must be finite and nonnegative.")
  group_counts <- table(metadata[[group_col]])
  # At least ten cells per fine group are required before probability inference.
  valid_groups <- names(group_counts)[group_counts >= min_cells]
  if (length(valid_groups) < 2) stop("fewer than two CellChat groups passed the minimum-cell filter")

  keep <- metadata[[group_col]] %in% valid_groups
  expr <- expr[, keep, drop = FALSE]
  metadata <- metadata[keep, , drop = FALSE]
  metadata[[group_col]] <- droplevels(factor(metadata[[group_col]]))

  cellchat <- createCellChat(object = expr, meta = metadata, group.by = group_col)
  cellchat@DB <- database
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(
    cellchat,
    type = "triMean",
    population.size = FALSE,
    nboot = 100L,
    seed.use = 1L,
    raw.use = TRUE
  )
  # raw.use selects unprojected signaling expression, not the input count scale.
  # triMean estimates robust group expression; cell abundance does not weight it.
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  cellchat
}


infer_cellchat_groups <- function(expr, metadata, diagnosis_col, database) {
  metadata <- align_cellchat_metadata(expr, metadata, c("fine", "origin", diagnosis_col))
  groups <- unique(data.table(origin = as.character(metadata$origin), diagnosis_group = as.character(metadata[[diagnosis_col]])))
  groups[, object_id := paste(origin, diagnosis_group, sep = "|")]
  if (anyDuplicated(groups$object_id)) stop("Group identifiers are not unique.")
  objects <- setNames(vector("list", nrow(groups)), groups$object_id)
  for (i in seq_len(nrow(groups))) {
    keep <- metadata$origin == groups$origin[i] & metadata[[diagnosis_col]] == groups$diagnosis_group[i]
    objects[[i]] <- run_cellchat_pipeline(expr[, keep, drop = FALSE], metadata[keep, , drop = FALSE], database)
  }
  list(objects = objects, groups = groups)
}


extract_spp1_communication <- function(inference) {
  rows <- lapply(seq_len(nrow(inference$groups)), function(i) {
    g <- inference$groups[i]
    obj <- inference$objects[[g$object_id]]
    if (is.null(obj)) stop("An inferred object is missing for ", g$object_id)
    prob <- obj@net$prob
    pval <- obj@net$pval
    if (is.null(prob) || is.null(pval) || !identical(dim(prob), dim(pval))) {
      stop("Inferred communication probability and p-value arrays are required.")
    }
    # A fully tested network can contain no significant links at alpha = 0.05.
    if (!any(is.finite(prob) & prob > 0 & is.finite(pval) & pval < 0.05)) {
      return(data.table(source = character(), target = character(), ligand = character(),
                        receptor = character(), prob = numeric(), pval = numeric(),
                        origin = character(), diagnosis_group = character()))
    }
    comm <- as.data.table(subsetCommunication(obj))
    comm <- comm[ligand == signaling_of_interest]
    comm[, `:=`(origin = g$origin, diagnosis_group = g$diagnosis_group)]
    comm
  })
  rbindlist(rows, fill = TRUE)
}


load_fine_colors <- function(colors, fine_levels) {
  if (length(colors) && (is.null(names(colors)) || anyDuplicated(names(colors)))) stop("Fine colors must be uniquely named.")
  colors <- colors[!is.na(colors) & nzchar(colors)]
  missing <- setdiff(fine_levels, names(colors))
  if (length(missing) > 0) {
    colors <- c(colors, setNames(grDevices::hcl.colors(length(missing), palette = "Dark 3"), missing))
  }
  colors[fine_levels]
}


build_spp1_role_matrix <- function(cellchat) {
  if (length(cellchat@netP$centr) == 0 || !signaling_of_interest %in% names(cellchat@netP$centr)) {
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  }
  centr <- cellchat@netP$centr[[signaling_of_interest]]
  if (is.null(centr)) stop("SPP1 centrality is unavailable in this inferred network.")
  role_mat <- matrix(unlist(centr), ncol = length(centr), byrow = FALSE)
  role_mat <- t(role_mat)
  rownames(role_mat) <- names(centr)
  colnames(role_mat) <- names(centr$outdeg)
  role_mat <- role_mat[c("outdeg", "indeg", "flowbet", "info"), , drop = FALSE]
  rownames(role_mat) <- c("Sender", "Receiver", "Mediator", "Influencer")
  # Within-role maximum scaling compares cell groups without mixing centrality units.
  row_max <- apply(role_mat, 1, max, na.rm = TRUE)
  row_max[!is.finite(row_max) | row_max == 0] <- 1
  sweep(role_mat, 1, row_max, "/", check.margin = FALSE)
}


prepare_spp1_lr_matrix <- function(cellchat, receiver_order) {
  # Each receptor column records the strongest detected sender probability per LR pair.
  comm <- as.data.table(subsetCommunication(cellchat, signaling = signaling_of_interest))
  if (nrow(comm) == 0) stop("no SPP1 communication rows were detected")
  if (!"interaction_name_2" %in% colnames(comm)) {
    comm[, interaction_name_2 := paste(ligand, receptor, sep = " - ")]
  }
  comm[, target := as.character(target)]
  comm[, interaction_name_2 := as.character(interaction_name_2)]
  lr_order <- unique(comm[order(interaction_name_2)]$interaction_name_2)
  lr_matrix <- matrix(
    NA_real_,
    nrow = length(lr_order),
    ncol = length(receiver_order),
    dimnames = list(lr_order, receiver_order)
  )
  summary <- comm[target %in% receiver_order, .(prob = max(prob, na.rm = TRUE)), by = .(interaction_name_2, target)]
  for (i in seq_len(nrow(summary))) {
    lr_matrix[summary$interaction_name_2[[i]], summary$target[[i]]] <- summary$prob[[i]]
  }
  lr_matrix
}


plot_spp1_combined_heatmap <- function(cellchat, fine_colors = character()) {
  receiver_order <- levels(cellchat@idents)
  if (lyve1_sender %in% receiver_order) receiver_order <- c(lyve1_sender, setdiff(receiver_order, lyve1_sender))

  role_mat <- build_spp1_role_matrix(cellchat)[, receiver_order, drop = FALSE]
  lr_mat <- prepare_spp1_lr_matrix(cellchat, receiver_order)
  fine_colors <- load_fine_colors(fine_colors, receiver_order)
  column_annotation <- HeatmapAnnotation(
    Fine = factor(receiver_order, levels = receiver_order),
    col = list(Fine = fine_colors),
    which = "column",
    show_annotation_name = FALSE,
    show_legend = FALSE,
    simple_anno_size = unit(2.0, "mm")
  )

  role_col_fun <- circlize::colorRamp2(c(0, 0.45, 1), c("#FFFFFF", "#31A354", "#005A32"))
  prob_range <- range(lr_mat, na.rm = TRUE)
  if (!all(is.finite(prob_range))) {
    prob_range <- c(0, 1)
  } else if (prob_range[1] == prob_range[2]) {
    prob_range <- c(0, max(1, prob_range[2]))
  }
  prob_col_fun <- circlize::colorRamp2(
    seq(prob_range[1], prob_range[2], length.out = 5),
    c("#08519C", "#6BAED6", "#FFFFFF", "#FB6A4A", "#B10026")
  )

  ht_role <- Heatmap(
    role_mat,
    name = "Importance",
    col = role_col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    height = unit(3.0, "cm"),
    column_title = "SPP1 signaling pathway network"
  )
  ht_lr <- Heatmap(
    lr_mat,
    name = "Commun.\nProb.",
    col = prob_col_fun,
    na_col = "white",
    bottom_annotation = column_annotation,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_labels = display_fine_label(receiver_order),
    column_names_rot = 45,
    height = unit(2.4, "cm")
  )
  list(plot = ht_role %v% ht_lr, role_matrix = role_mat, lr_matrix = lr_mat)
}


plot_spp1_pair_fraction_bubble <- function(spp1_table) {
  spp1 <- copy(as.data.table(spp1_table))
  if ("diagnosis_group" %in% colnames(spp1) && !"b_diagnosis_new" %in% colnames(spp1)) {
    setnames(spp1, "diagnosis_group", "b_diagnosis_new")
  }
  required <- c("source", "target", "ligand", "receptor", "prob", "pval", "origin", "b_diagnosis_new")
  missing <- setdiff(required, colnames(spp1))
  if (length(missing) > 0) stop("SPP1 table missing columns: ", paste(missing, collapse = ", "))
  spp1 <- spp1[ligand == signaling_of_interest]
  if (nrow(spp1) == 0) stop("no SPP1 rows are available for pair-fraction bubble plotting")
  receptor_levels <- sort(unique(as.character(spp1$receptor)))
  origin_levels <- c("CSF", "Blood")
  origin_levels <- c(origin_levels[origin_levels %in% unique(spp1$origin)], setdiff(sort(unique(spp1$origin)), origin_levels))
  diagnosis_levels <- c("CN_A-", "Preclinical_AD", "MCI_due_to_AD", "Dementia_due_to_AD", "Non-AD")
  diagnosis_levels <- c(diagnosis_levels[diagnosis_levels %in% unique(spp1$b_diagnosis_new)], setdiff(sort(unique(spp1$b_diagnosis_new)), diagnosis_levels))

  # The published denominator is conditional on SPP1-active source/target groups.
  # It is not the complete space of all cell-group pairs retained for inference.
  pair_space <- spp1[, .(
    n_possible_sources = uniqueN(source),
    n_possible_targets = uniqueN(target),
    possible_source_target_rows = uniqueN(source) * uniqueN(target)
  ), by = .(origin, b_diagnosis_new)]
  observed <- spp1[, .(
    prob = max(prob, na.rm = TRUE),
    detected_source_target_rows = uniqueN(paste(source, target, sep = "\r")),
    raw_rows = .N
  ), by = .(origin, b_diagnosis_new, receptor)]
  grid_dt <- CJ(origin = origin_levels, b_diagnosis_new = diagnosis_levels, receptor = receptor_levels)
  plot_dt <- merge(grid_dt, pair_space, by = c("origin", "b_diagnosis_new"), all.x = TRUE)
  plot_dt <- merge(plot_dt, observed, by = c("origin", "b_diagnosis_new", "receptor"), all.x = TRUE)
  plot_dt[, detected := !is.na(detected_source_target_rows)]
  plot_dt[is.na(detected_source_target_rows), detected_source_target_rows := 0L]
  plot_dt[is.na(raw_rows), raw_rows := 0L]
  plot_dt[is.na(possible_source_target_rows), possible_source_target_rows := 0L]
  plot_dt[, detected_pair_fraction := fifelse(detected & possible_source_target_rows > 0, detected_source_target_rows / possible_source_target_rows, NA_real_)]
  plot_dt[, `:=`(
    origin = factor(origin, levels = origin_levels),
    b_diagnosis_new = factor(b_diagnosis_new, levels = diagnosis_levels),
    receptor = factor(receptor, levels = rev(receptor_levels))
  )]

  prob_range <- range(plot_dt$prob[plot_dt$detected], na.rm = TRUE)
  if (!any(plot_dt$detected)) stop("no detected SPP1 origin-diagnosis-receptor combinations are available")
  if (!all(is.finite(prob_range)) || prob_range[1] == prob_range[2]) prob_range <- c(0, max(1, plot_dt$prob, na.rm = TRUE))

  p <- ggplot(plot_dt, aes(x = b_diagnosis_new, y = receptor)) +
    geom_point(data = plot_dt[detected == FALSE], shape = 21, size = 2.2, color = "grey75", fill = "white", stroke = 0.35) +
    geom_point(data = plot_dt[detected == TRUE], aes(size = detected_pair_fraction, color = prob), alpha = 0.96) +
    scale_color_gradientn(name = "Communication\nprobability", colors = c("#F7FBFF", "#D8ECF7", "#FB6A4A", "#B10026"), limits = prob_range) +
    scale_size_continuous(name = "Detected pair fraction\n(SPP1-active groups)", range = c(2.0, 7.4), limits = c(0, max(plot_dt$detected_pair_fraction, na.rm = TRUE))) +
    facet_grid(. ~ origin) +
    theme_bw(base_family = "Arial", base_size = 8.5) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      strip.background = element_blank(),
      strip.text = element_text(size = 9),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.2)
    )
  list(plot = p, data = plot_dt)
}


extract_lyve1_network <- function(inference, sender = lyve1_sender) {
  receiver_rows <- list()
  pathway_rows <- list()
  for (i in seq_len(nrow(inference$groups))) {
    g <- inference$groups[i]
    obj <- inference$objects[[g$object_id]]
    if (is.null(obj)) stop("An inferred object is missing for ", g$object_id)
    if (!sender %in% levels(obj@idents)) next
    pathways <- as.character(obj@netP$pathways)
    for (j in seq_along(pathways)) {
      values <- obj@netP$prob[sender, , j]
      pathway_rows[[length(pathway_rows) + 1L]] <- data.table(
        origin = g$origin, diagnosis_group = g$diagnosis_group,
        feature = pathways[j], label = pathways[j],
        probability = sum(values), n_receiver_positive = sum(values > 0)
      )
      if (pathways[j] == signaling_of_interest) {
        receiver_rows[[length(receiver_rows) + 1L]] <- data.table(
          origin = g$origin, diagnosis_group = g$diagnosis_group,
          feature = names(values), label = display_fine_label(names(values)),
          probability = as.numeric(values)
        )
      }
    }
  }
  empty <- data.table(origin = character(), diagnosis_group = character(),
                      feature = character(), label = character(), probability = numeric())
  list(
    receiver = if (length(receiver_rows)) rbindlist(receiver_rows) else copy(empty),
    pathway = if (length(pathway_rows)) rbindlist(pathway_rows) else copy(empty)
  )
}


plot_lyve1_chord <- function(network, origin_name, mode = c("receiver", "pathway"),
                             top_n = NULL, fine_colors = character()) {
  mode <- match.arg(mode)
  if (is.null(top_n)) top_n <- if (mode == "receiver") 20L else 30L
  if (length(top_n) != 1L || !is.finite(top_n) || top_n < 1 || top_n != as.integer(top_n)) {
    stop("top_n must be a positive integer.")
  }
  dt <- copy(network[[mode]][origin == origin_name & probability > 0])
  if (!nrow(dt)) stop("No positive LYVE1 ", mode, " links for ", origin_name)
  # Ranking uses the largest group-specific probability; pathway views retain SPP1.
  score <- dt[, .(maximum = max(probability), total = sum(probability)), by = feature]
  setorder(score, -maximum, -total)
  keep <- head(score$feature, top_n)
  if (mode == "pathway") keep <- unique(c(keep, signaling_of_interest))
  dt <- dt[feature %in% keep]
  dt[, value := probability * 1e6]

  groups <- unique(dt$diagnosis_group)
  group_colors <- c("CN_A-" = "#7AA6C2", "Biological_AD" = "#CB181D", "Non-AD" = "#8F8F8F")
  group_colors <- load_fine_colors(group_colors, groups)
  targets <- dt[, .(total = sum(value), maximum = max(value)), by = .(feature, label)]
  if (mode == "receiver") setorder(targets, -total) else setorder(targets, maximum)
  source_ids <- paste0("source:", groups)
  target_ids <- paste0("target:", targets$feature)
  source_map <- setNames(source_ids, groups)
  target_map <- setNames(target_ids, targets$feature)
  sector_order <- c(source_ids, target_ids)
  source_labels <- gsub("_", " ", groups)
  if (mode == "receiver") source_labels <- paste(source_labels, "LYVE1", sep = "\n")
  sector_labels <- setNames(c(source_labels, targets$label), sector_order)
  target_colors <- if (mode == "receiver") {
    load_fine_colors(fine_colors, targets$feature)
  } else {
    setNames(ifelse(targets$feature == signaling_of_interest, "#FEE8C8", "#EDEDED"), targets$feature)
  }
  sector_colors <- setNames(c(unname(group_colors), unname(target_colors)), sector_order)
  source_totals <- dt[, .(total = sum(value)), by = diagnosis_group]
  widths <- c(setNames(source_totals$total, source_map[source_totals$diagnosis_group]),
              setNames(targets$total, target_ids))
  limits <- range(dt$value)
  if (diff(limits) == 0) limits <- c(0, max(1, limits[2]))
  probability_colors <- circlize::colorRamp2(
    c(limits[1], mean(limits), limits[2]), c("#FEE8C8", "#F46D43", "#7F0000")
  )

  op <- par(no.readonly = TRUE)
  on.exit({ circlize::circos.clear(); par(op) }, add = TRUE)
  par(mar = c(1.2, 1.2, 3.0, 1.2), family = "Arial")
  circos.clear()
  circos.par(start.degree = 102,
             cell.padding = c(0.02, 0, 0.02, 0),
             canvas.xlim = c(-1.2, 1.2), canvas.ylim = c(-1.15, 1.15),
             gap.after = c(rep(3.5, length(source_ids) - 1), 16,
                           rep(2.0, length(target_ids) - 1), 22))
  circos.initialize(sector_order, xlim = cbind(0, pmax(widths[sector_order], 1)))
  circos.trackPlotRegion(
    ylim = c(0, 1.6), track.height = 0.11, bg.border = NA,
    panel.fun = function(x, y) {
      sector <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      highlight <- mode == "pathway" && sector == paste0("target:", signaling_of_interest)
      circos.rect(xlim[1], 0, xlim[2], 1, col = sector_colors[[sector]],
                  border = if (highlight) "black" else darken_hex(sector_colors[[sector]], 0.56),
                  lwd = if (highlight) 2.2 else 1.2)
      circos.text(mean(xlim), 1.25, sector_labels[[sector]], facing = "clockwise",
                  niceFacing = TRUE, adj = c(0, 0.5), cex = 0.44, family = "Arial")
    }
  )
  offsets <- setNames(rep(0, length(sector_order)), sector_order)
  setorder(dt, diagnosis_group, -value)
  for (i in seq_len(nrow(dt))) {
    src <- source_map[[dt$diagnosis_group[i]]]
    tgt <- target_map[[dt$feature[i]]]
    value <- dt$value[i]
    src_range <- c(offsets[[src]], offsets[[src]] + value)
    tgt_range <- c(offsets[[tgt]], offsets[[tgt]] + value)
    offsets[[src]] <- src_range[2]
    offsets[[tgt]] <- tgt_range[2]
    link_color <- if (mode == "receiver") probability_colors(value) else group_colors[[dt$diagnosis_group[i]]]
    circos.link(src, src_range, tgt, tgt_range,
                col = adjustcolor(link_color, alpha.f = 0.82),
                border = NA, directional = 1, arr.type = "big.arrow", arr.length = 0.018)
  }
  grid.text(
    paste0(origin_name, ": ncMono_05_LYVE1 ", if (mode == "receiver") "SPP1 receivers" else "outgoing signals"),
    x = 0.5, y = 0.965, gp = gpar(fontfamily = "Arial", fontsize = 10.5, fontface = "bold")
  )
  invisible(dt)
}


select_cellchat_sample <- function(metadata, sampled_cells) {
  # The biological-diagnosis analysis uses floor(0.30 * n) cells per orig.ident/fine stratum.
  # The frozen barcode selection preserves the original NumPy PCG64 seed 20260525 draw.
  if (!all(c("orig.ident", "fine") %in% names(metadata))) stop("Sampling requires orig.ident and fine.")
  if (!length(sampled_cells) || anyNA(sampled_cells) || anyDuplicated(sampled_cells) ||
      !all(sampled_cells %in% rownames(metadata))) stop("Sampled barcodes must be a unique subset of metadata.")
  strata <- data.table(metadata[, c("orig.ident", "fine"), drop = FALSE])
  strata[, selected := rownames(metadata) %in% sampled_cells]
  sizes <- strata[, .(n_full = .N, n_sample = sum(selected)), by = .(orig.ident, fine)]
  if (any(sizes$n_sample != floor(0.30 * sizes$n_full))) {
    stop("Sampled barcodes do not match 30% sampling within orig.ident/fine strata.")
  }
  which(rownames(metadata) %in% sampled_cells)
}


run_cellchat_analysis <- function(expr, metadata, sampled_cells, expression_scale,
                                   database = CellChatDB.human, fine_colors = character()) {
  # Input is log-normalized RNA expression, not counts or integrated residuals.
  # https://github.com/jinworks/CellChat/blob/main/man/createCellChat.Rd
  expression_scale <- match.arg(expression_scale, "lognormalized")
  metadata <- align_cellchat_metadata(
    expr, metadata, c("fine", "origin", "b_diagnosis_new", "Biological_diagnosis", "orig.ident")
  )
  values <- if (inherits(expr, "sparseMatrix")) expr@x else as.vector(expr)
  if (any(!is.finite(values)) || any(values < 0)) stop("Expression must be finite and nonnegative.")
  # The archived analyses use the complete human database, including all signaling classes.
  # Biological_diagnosis is the curated CN_A-/Biological_AD/Non-AD annotation.
  biological_levels <- c("CN_A-", "Biological_AD", "Non-AD")
  if (!all(metadata$Biological_diagnosis %in% biological_levels)) stop("Unexpected Biological_diagnosis labels.")
  sample_index <- select_cellchat_sample(metadata, sampled_cells)
  stage <- infer_cellchat_groups(expr, metadata, "b_diagnosis_new", database)
  biological <- infer_cellchat_groups(
    expr[, sample_index, drop = FALSE], metadata[sample_index, , drop = FALSE],
    "Biological_diagnosis", database
  )
  spp1 <- extract_spp1_communication(stage)
  heatmap_id <- stage$groups[origin == "CSF" & diagnosis_group == "Preclinical_AD", object_id]
  if (length(heatmap_id) != 1L) stop("The CSF Preclinical_AD group is required for the target heatmap.")
  list(
    stage = stage, biological = biological, spp1 = spp1,
    expression_scale = expression_scale,
    sampled_cells = colnames(expr)[sample_index],
    heatmap = plot_spp1_combined_heatmap(stage$objects[[heatmap_id]], fine_colors),
    bubble = plot_spp1_pair_fraction_bubble(spp1),
    lyve1_network = extract_lyve1_network(biological)
  )
}
