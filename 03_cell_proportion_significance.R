# Significance testing for single-cell proportions
#
# This script contains no file I/O or plotting. Source it, then pass a metadata
# data frame to run_proportion_significance().
#
# Final comparisons:
# - AD Continuum: three adjacent b_diagnosis contrasts.
# - Dementia Course: two adjacent c_diagnosis contrasts.

library(dplyr)
library(tidyr)
library(stringr)
library(emmeans)

AD_PROGRESSION_LEVELS <- c(
  "CN A-",
  "Preclinical AD",
  "MCI due to AD",
  "Dementia due to AD"
)

CLINICAL_STAGE_LEVELS <- c("CN", "MCI", "Dementia")

require_columns <- function(data, columns) {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

prepare_significance_metadata <- function(
  metadata) {
  require_columns(
    metadata,
    c(
      "patientid", "orig.ident", "origin", "celltype", "fine",
      "b_diagnosis", "c_diagnosis", "age", "sex"
    )
  )

  metadata %>%
    transmute(
      patientid,
      sampleid = orig.ident,
      sample_type = str_trim(as.character(origin)),
      celltype = as.character(celltype),
      fine = str_remove_all(str_trim(as.character(fine)), "\\s+"),
      DX = as.character(b_diagnosis),
      b_diagnosis = as.character(b_diagnosis),
      c_diagnosis = str_trim(as.character(c_diagnosis)),
      age = as.numeric(age),
      sex = factor(sex)
    ) %>%
    filter(
      !is.na(sampleid), !is.na(sample_type), !is.na(celltype), !is.na(fine),
      !is.na(DX), !is.na(age), !is.na(sex)
    )
}

default_comparisons <- function() {
  list(
    list(
      category = "1. AD Progression",
      category_short = "AD_Progression",
      name = "AD_trajectory_adjacent_pairs",
      groups = AD_PROGRESSION_LEVELS,
      contrast_type = "biological_stage",
      requested_pairs = list(
        c("Preclinical AD", "CN A-"),
        c("MCI due to AD", "Preclinical AD"),
        c("Dementia due to AD", "MCI due to AD")
      ),
      desc = "AD trajectory adjacent comparisons",
      priority = "high"
    ),
    list(
      category = "2. Clinical Stage",
      category_short = "Clinical_Stage",
      name = "stage_adjacent_pairs",
      groups = CLINICAL_STAGE_LEVELS,
      contrast_type = "clinical_stage",
      requested_pairs = list(c("MCI", "CN"), c("Dementia", "MCI")),
      desc = "Clinical stage adjacent comparisons (c_diagnosis)",
      priority = "medium"
    )
  )
}

aggregate_proportion_counts <- function(data, analysis_level, target) {
  if (analysis_level == "celltype") {
    sample_data <- data %>%
      distinct(sampleid, DX, c_diagnosis, age, sex, b_diagnosis)
    denominators <- data %>%
      count(sampleid, DX, name = "total_cells")
    numerators <- data %>%
      filter(celltype == target) %>%
      count(sampleid, DX, name = "n_cells")
    parent_celltype <- NA_character_
  } else {
    parent_celltypes <- data %>%
      filter(fine == target) %>%
      distinct(celltype) %>%
      pull(celltype)

    if (length(parent_celltypes) != 1) {
      stop("Fine subtype must map to exactly one parent cell type.", call. = FALSE)
    }

    parent_celltype <- parent_celltypes[[1]]
    parent_data <- filter(data, celltype == parent_celltype)
    sample_data <- parent_data %>%
      distinct(sampleid, DX, c_diagnosis, age, sex, b_diagnosis)
    denominators <- parent_data %>%
      count(sampleid, DX, name = "total_cells")
    numerators <- parent_data %>%
      filter(fine == target) %>%
      count(sampleid, DX, name = "n_cells")
  }

  counts <- sample_data %>%
    left_join(denominators, by = c("sampleid", "DX")) %>%
    left_join(numerators, by = c("sampleid", "DX")) %>%
    mutate(n_cells = replace_na(n_cells, 0L), sex = factor(sex))

  list(data = counts, parent_celltype = parent_celltype)
}

prepare_comparison <- function(data, comparison) {
  type <- comparison$contrast_type

  if (type == "clinical_stage") {
    prepared <- data %>%
      mutate(
        comparison_group = if_else(
          c_diagnosis %in% c("CN", "MCI", "Dementia"),
          c_diagnosis,
          NA_character_
        ),
        comparison_group = factor(
          comparison_group,
          levels = c("CN", "MCI", "Dementia")
        )
      ) %>%
      filter(!is.na(comparison_group))
  } else if (type == "biological_stage") {
    prepared <- data %>%
      filter(DX %in% comparison$groups) %>%
      mutate(comparison_group = factor(DX, levels = comparison$groups))
  } else {
    stop("Unknown contrast_type: ", type, call. = FALSE)
  }

  prepared
}

parse_contrast <- function(contrast) {
  clean <- str_remove_all(as.character(contrast), "[()]")
  parts <- str_split_fixed(clean, " - ", 2)
  tibble(
    group1 = str_trim(parts[, 1]),
    group2 = str_trim(parts[, 2])
  )
}

summarise_groups <- function(data) {
  data %>%
    mutate(proportion = n_cells / total_cells) %>%
    group_by(comparison_group) %>%
    summarise(
      n_samples = n(),
      mean_proportion = mean(proportion),
      sd_proportion = sd(proportion),
      median_proportion = median(proportion),
      n_zero_samples = sum(n_cells == 0),
      pct_zero_samples = 100 * n_zero_samples / n_samples,
      mean_n_cells = mean(n_cells),
      mean_total_cells = mean(total_cells),
      .groups = "drop"
    ) %>%
    mutate(comparison_group = as.character(comparison_group))
}

summarise_composition <- function(data) {
  data %>%
    count(comparison_group, b_diagnosis, name = "n") %>%
    group_by(comparison_group) %>%
    summarise(
      composition = paste0(b_diagnosis, "(", n, ")", collapse = "; "),
      .groups = "drop"
    ) %>%
    mutate(comparison_group = as.character(comparison_group))
}

add_group_details <- function(result, data) {
  contrast_groups <- parse_contrast(result$contrast)
  result <- bind_cols(result, contrast_groups)
  summaries <- summarise_groups(data)
  compositions <- summarise_composition(data)

  group_details <- summaries %>%
    left_join(compositions, by = "comparison_group")

  result %>%
    left_join(group_details, by = c("group1" = "comparison_group")) %>%
    rename_with(~ paste0("group1_", .x), all_of(names(group_details)[-1])) %>%
    left_join(group_details, by = c("group2" = "comparison_group")) %>%
    rename_with(~ paste0("group2_", .x), all_of(names(group_details)[-1]))
}

keep_requested_pairs <- function(result, requested_pairs) {
  if (is.null(requested_pairs)) return(result)
  requested_keys <- vapply(
    requested_pairs,
    function(x) paste(x[[1]], x[[2]], sep = "\r"),
    character(1)
  )
  result %>%
    filter(paste(group1, group2, sep = "\r") %in% requested_keys)
}

fit_one_comparison <- function(
    counts,
    comparison,
    tissue,
    analysis_level,
    target,
    parent_celltype,
    min_samples_per_group) {
  data <- prepare_comparison(counts, comparison)
  sample_counts <- table(data$comparison_group)
  sample_counts <- sample_counts[sample_counts > 0]

  if (length(sample_counts) < 2 || any(sample_counts < min_samples_per_group)) {
    return(list(
      result = NULL,
      reason = "Insufficient samples per group",
      details = paste(names(sample_counts), sample_counts, sep = "=", collapse = ", ")
    ))
  }

  model <- glm(
    cbind(n_cells, total_cells - n_cells) ~ comparison_group + age + sex,
    data = data,
    family = quasibinomial()
  )
  emmeans_result <- emmeans(model, ~ comparison_group)
  # Use unadjusted planned contrasts here. BH correction is applied once,
  # after retaining only the prespecified adjacent comparisons.
  result <- as.data.frame(
    pairs(emmeans_result, reverse = TRUE, adjust = "none")
  )

  if ("z.ratio" %in% names(result)) {
    result$statistic <- result$z.ratio
  } else if ("t.ratio" %in% names(result)) {
    result$statistic <- result$t.ratio
  } else {
    result$statistic <- NA_real_
  }

  result <- result %>%
    select(-any_of(c("z.ratio", "t.ratio"))) %>%
    add_group_details(data)
  result <- keep_requested_pairs(result, comparison$requested_pairs)
  result <- result %>%
    mutate(
      category = comparison$category,
      category_short = comparison$category_short,
      comparison_name = comparison$name,
      comparison_desc = comparison$desc,
      priority = comparison$priority,
      tissue = tissue,
      level = analysis_level,
      celltype = if_else(analysis_level == "celltype", target, NA_character_),
      parent_celltype = if_else(
        analysis_level == "fine",
        parent_celltype,
        NA_character_
      ),
      fine_subtype = if_else(analysis_level == "fine", target, NA_character_),
      covariates = "age + sex",
      dispersion = summary(model)$dispersion,
      n_samples_total = nrow(data),
      proportion_definition = if (analysis_level == "celltype") {
        paste0(target, " / all cells in each sample")
      } else {
        paste0(target, " / all ", parent_celltype, " cells in each sample")
      }
    ) %>%
    relocate(
      category, category_short, comparison_name, comparison_desc, priority,
      tissue, level, celltype, parent_celltype, fine_subtype,
      contrast, group1, group2, estimate, SE, df, statistic, p.value
    )

  list(result = result, reason = NULL, details = NULL)
}

apply_fdr <- function(results) {
  results %>%
    mutate(
      fdr_family = case_when(
        category == "1. AD Progression" ~ "AD Continuum",
        category == "2. Clinical Stage" ~ "Dementia Course",
        TRUE ~ NA_character_
      )
    ) %>%
    group_by(category, tissue) %>%
    mutate(fdr_category_tissue = p.adjust(p.value, method = "BH")) %>%
    ungroup() %>%
    mutate(
      # Backward-compatible alias for older plotting code.
      fdr_family_tissue = fdr_category_tissue,
      sig_level = case_when(
        fdr_category_tissue < 0.001 ~ "***",
        fdr_category_tissue < 0.01 ~ "**",
        fdr_category_tissue < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      significant_raw = p.value < 0.05,
      significant_category_tissue = fdr_category_tissue < 0.05,
      significant_family_tissue = significant_category_tissue,
      contrast_vs = str_replace_all(str_remove_all(contrast, "[()]"), " - ", " vs ")
    )
}

# The returned significant table uses BH-FDR < 0.05 within each final
# analysis category and tissue compartment.
run_proportion_significance <- function(
    metadata,
    tissues = c("CSF", "Blood"),
    analysis_levels = c("celltype", "fine"),
    comparisons = default_comparisons(),
    min_samples_per_group = 3) {
  clean <- prepare_significance_metadata(metadata)
  results <- list()
  skipped <- list()
  result_index <- 0L
  skip_index <- 0L

  for (tissue in tissues) {
    tissue_data <- filter(clean, sample_type == tissue)

    for (analysis_level in analysis_levels) {
      targets <- if (analysis_level == "celltype") {
        unique(tissue_data$celltype)
      } else {
        unique(tissue_data$fine)
      }

      for (target in targets) {
        aggregated <- tryCatch(
          aggregate_proportion_counts(tissue_data, analysis_level, target),
          error = identity
        )

        if (inherits(aggregated, "error")) {
          skip_index <- skip_index + 1L
          skipped[[skip_index]] <- tibble(
            tissue = tissue,
            level = analysis_level,
            target = target,
            comparison = NA_character_,
            reason = "Aggregation error",
            details = conditionMessage(aggregated)
          )
          next
        }

        for (comparison in comparisons) {
          fitted <- tryCatch(
            fit_one_comparison(
              aggregated$data,
              comparison,
              tissue,
              analysis_level,
              target,
              aggregated$parent_celltype,
              min_samples_per_group
            ),
            error = function(error) {
              list(
                result = NULL,
                reason = "Model error",
                details = conditionMessage(error)
              )
            }
          )

          if (is.null(fitted$result)) {
            skip_index <- skip_index + 1L
            skipped[[skip_index]] <- tibble(
              tissue = tissue,
              level = analysis_level,
              target = target,
              comparison = comparison$name,
              reason = fitted$reason,
              details = fitted$details
            )
          } else {
            result_index <- result_index + 1L
            results[[result_index]] <- fitted$result
          }
        }
      }
    }
  }

  if (length(results) == 0) {
    stop("No models were successfully fitted.", call. = FALSE)
  }

  all_results <- apply_fdr(bind_rows(results))
  skipped_results <- if (length(skipped) == 0) tibble() else bind_rows(skipped)

  list(
    results = all_results,
    significant = filter(all_results, significant_category_tissue),
    skipped = skipped_results
  )
}
