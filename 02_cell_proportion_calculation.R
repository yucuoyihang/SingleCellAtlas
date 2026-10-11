library(dplyr)
library(stringr)

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

prepare_proportion_metadata <- function(metadata) {
  require_columns(
    metadata,
    c("patientid", "origin", "celltype", "fine", "b_diagnosis", "c_diagnosis")
  )

  metadata %>%
    mutate(
      fine = str_remove_all(str_trim(fine), "\\s+"),
      DX = as.character(b_diagnosis),
      # Clinical stage is defined independently of biological diagnosis.
      DX_stage = str_trim(as.character(c_diagnosis)),
      sample_type = str_trim(as.character(origin)),
      celltype = as.character(celltype)
    ) %>%
    filter(
      !is.na(patientid), !is.na(DX), !is.na(celltype), !is.na(fine),
      !is.na(sample_type), sample_type != ""
    )
}

# Major-cell proportion uses all cells in each sample_type/diagnosis group.
# Fine-subtype proportion uses all cells in its parent major-cell group.
calculate_proportions_both_levels <- function(data, diagnosis_col = "DX") {
  require_columns(data, c("sample_type", diagnosis_col, "celltype", "fine"))

  celltype <- data %>%
    filter(!is.na(.data[[diagnosis_col]])) %>%
    transmute(
      sample_type,
      diagnosis = .data[[diagnosis_col]],
      celltype
    ) %>%
    count(sample_type, diagnosis, celltype, name = "n_celltype") %>%
    group_by(sample_type, diagnosis) %>%
    mutate(
      n_total = sum(n_celltype),
      proportion = 100 * n_celltype / n_total
    ) %>%
    ungroup() %>%
    mutate(group_id = paste(sample_type, diagnosis, sep = "_")) %>%
    rename(DX = diagnosis)

  fine <- data %>%
    filter(!is.na(.data[[diagnosis_col]])) %>%
    transmute(
      sample_type,
      diagnosis = .data[[diagnosis_col]],
      celltype,
      fine
    ) %>%
    count(sample_type, diagnosis, celltype, fine, name = "n_fine") %>%
    group_by(sample_type, diagnosis, celltype) %>%
    mutate(
      n_celltype = sum(n_fine),
      proportion = 100 * n_fine / n_celltype
    ) %>%
    ungroup() %>%
    mutate(group_id = paste(sample_type, diagnosis, sep = "_")) %>%
    rename(DX = diagnosis)

  list(celltype = celltype, fine = fine)
}

find_paired_patients <- function(data, sample_types = c("CSF", "Blood")) {
  require_columns(data, c("patientid", "sample_type"))

  data %>%
    distinct(patientid, sample_type) %>%
    filter(sample_type %in% sample_types) %>%
    count(patientid, name = "n_sample_types") %>%
    filter(n_sample_types == length(sample_types)) %>%
    pull(patientid)
}

# Returns the four proportion sets used in the original analysis plus the IDs
# included in the paired CSF/Blood subset.
calculate_all_proportions <- function(metadata) {
  clean <- prepare_proportion_metadata(metadata)
  paired_patient_ids <- find_paired_patients(clean)
  paired <- filter(clean, patientid %in% paired_patient_ids)

  list(
    progression_all = calculate_proportions_both_levels(
      filter(clean, DX %in% AD_PROGRESSION_LEVELS),
      "DX"
    ),
    stage_all = calculate_proportions_both_levels(
      filter(clean, DX_stage %in% CLINICAL_STAGE_LEVELS),
      "DX_stage"
    ),
    progression_paired = calculate_proportions_both_levels(
      filter(paired, DX %in% AD_PROGRESSION_LEVELS),
      "DX"
    ),
    stage_paired = calculate_proportions_both_levels(
      filter(paired, DX_stage %in% CLINICAL_STAGE_LEVELS),
      "DX_stage"
    ),
    paired_patient_ids = paired_patient_ids
  )
}
