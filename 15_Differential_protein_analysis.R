# Differential protein analysis ================================================

library(data.table)
library(dplyr)

# Load QCed proteomic data, protein annotation, and clinical metadata
load("pro_QCed.RData")
load("pro_ANNO.RData") # information for all proteins
clinical <- fread("clinical_metadata.tsv")

covariates <- c("Age", "Sex")

# Select CSF or serum proteomic data
tissue <- "CSF"  # "CSF" or "serum"
pro_final <- get(paste0("pro_", tissue, "_qc")) %>%
  merge(clinical, by = "IDENTIFIER_1", all.x = TRUE) %>% 
  mutate(across(all_of(pro_ANNO$Target_run), scale))

# Define pairwise comparisons
Clinical_pairs <- list(
  c("MCI", "CN"),
  c("Dementia", "CN"),
  c("Dementia", "MCI")
)

Biological_pairs <- list(
  c("Preclinical_AD", "CN_A-"),
  c("MCI_due_to_AD", "CN_A-"),
  c("Dementia_due_to_AD", "CN_A-"),
  c("MCI_due_to_AD", "Preclinical_AD"),
  c("Dementia_due_to_AD", "Preclinical_AD"),
  c("Dementia_due_to_AD", "MCI_due_to_AD"),
  c("Non-AD", "CN_A-")
)

BiologicalAD_pairs <- list(
  c("Biological_AD", "CN_A-"),
  c("Biological_AD", "Non-AD")
)

# Example comparison
diag <- "BiologicalAD"
pair <- BiologicalAD_pairs[[1]]
case_group <- pair[1]
control_group <- pair[2]

pro_final_temp <- pro_final %>%
  filter(.data[[paste0(diag, "_diagnosis")]] %in% pair) %>%
  mutate(group = ifelse(.data[[paste0(diag, "_diagnosis")]] == case_group, 1, 0))

# Age- and sex-adjusted logistic regression for each aptamer
result_binary <- lapply(pro_ANNO$Target_run, function(n) {

  FML <- as.formula(paste(
    "group ~", n,
    if(length(covariates)) paste("+", paste(covariates, collapse = " + ")) else ""
  ))

  fit <- glm(FML, data = pro_final_temp, family = binomial)
  sfit <- summary(fit)

  data.frame(
    case = case_group,
    control = control_group,
    case_No = sum(fit$y == 1, na.rm = TRUE),
    control_No = sum(fit$y == 0, na.rm = TRUE),
    Target_run = n,
    BETA = sfit$coefficients[2, 1],
    SE = sfit$coefficients[2, 2],
    z = sfit$coefficients[2, 3],
    Pvalue = sfit$coefficients[2, 4]
  )
})

result_final <- bind_rows(result_binary) %>%
  merge(pro_ANNO %>% select(Target_run, EntrezGeneSymbol, Target, TargetFullName),
        by = "Target_run", all.x = TRUE) %>%
  arrange(Pvalue) %>%
  select(case, control, case_No, control_No, Target_run, EntrezGeneSymbol,
         Target, TargetFullName, everything())

fwrite(result_final, paste0(tissue, "_", paste(pair, collapse = "-vs-"), "_logistic.tsv"),
       sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")