# Calculating DEG

library(Seurat)
library(MAST)
library(tidyverse)

comparison_list <- list(
  list(name = "Biological_AD vs CN_A-", case = c("Dementia_due_to_AD", "MCI_due_to_AD", "Preclinical_AD"), control = c("CN_A-")),
  list(name = "Preclinical_AD vs CN_A-", case = c("Preclinical_AD"), control = c("CN_A-")),
  list(name = "MCI_due_to_AD vs CN_A-", case = c("MCI_due_to_AD"), control = c("CN_A-")),
  list(name = "MCI_due_to_AD vs Preclinical_AD", case = c("MCI_due_to_AD"), control = c("Preclinical_AD")),
  list(name = "Dementia_due_to_AD vs Preclinical_AD", case = c("Dementia_due_to_AD"), control = c("Preclinical_AD")),
  list(name = "Dementia_due_to_AD vs MCI_due_to_AD", case = c("Dementia_due_to_AD"), control = c("MCI_due_to_AD")),
  list(name = "Non_AD vs CN_A-", case = c("Dementia_due_to_non-AD", "MCI_due_to_non-AD"), control = c("CN_A-")),
  list(name = "Biological_AD vs Non-AD", case = c("Dementia_due_to_AD", "MCI_due_to_AD", "Preclinical_AD"), control = c("Dementia_due_to_non-AD", "MCI_due_to_non-AD")),
  list(name = "Dementia vs MCI", case = c("Dementia_due_to_AD", "Dementia_due_to_non-AD"), control = c("MCI_due_to_AD", "MCI_due_to_non-AD")),
  list(name = "Dementia vs CN", case = c("Dementia_due_to_AD", "Dementia_due_to_non-AD"), control = c("CN_A-")),
  list(name = "MCI vs CN", case = c("MCI_due_to_AD", "MCI_due_to_non-AD", "Dementia_due_to_non-AD"), control = c("CN_A-"))
)


run_deg_pipeline <- function(seurat_full, target_origin = "Blood", min_cells = 10, min_patients = 3) {
  
  seurat_origin <- subset(seurat_full, subset = origin == target_origin)
  
  for (j in seq_along(comparison_list)) {
    comp <- comparison_list[[j]]
    cat(sprintf("\n[Group %d/%d] Running: %s\n", j, length(comparison_list), comp$name))
    
    all_diagnoses <- c(comp$case, comp$control)
    seurat_sub <- subset(seurat_origin, subset = b_diagnosis %in% all_diagnoses)
    
    if (ncol(seurat_sub) == 0) {
      cat("  [Warning] Low number of cells in the comparision group and passing \n")
      next
    }
    
    seurat_sub$group_compare <- ifelse(seurat_sub$b_diagnosis %in% comp$case, "Case", "Control")
    seurat_sub$age <- as.numeric(seurat_sub$age)
    
    cell_patient_qc <- seurat_sub@meta.data %>%
      group_by(fine, group_compare) %>%
      summarise(
        n_cells = n(),
        n_patients = n_distinct(orig.ident),
        .groups = "drop"
      ) %>%
      pivot_wider(
        names_from = group_compare, 
        values_from = c(n_cells, n_patients), 
        values_fill = 0
      ) %>%
      filter(
        n_cells_Case >= min_cells & n_cells_Control >= min_cells,
        n_patients_Case >= min_patients & n_patients_Control >= min_patients
      )
    
    valid_celltypes <- cell_patient_qc$fine
    
    if (length(valid_celltypes) == 0) {
      cat("  [Pass] No cell types with enough counts\n")
      next
    }
    
    deg_results_list <- map(seq_along(valid_celltypes), function(i) {
      ct <- valid_celltypes[i]
      cat(sprintf("  -> Analyzing [%d/%d]: %s\n", i, length(valid_celltypes), ct))
      
      seurat_ct <- subset(seurat_sub, subset = fine == ct)
      
      tryCatch({
        markers <- FindMarkers(
          seurat_ct,
          ident.1 = "Case",
          ident.2 = "Control",
          group.by = "group_compare",
          logfc.threshold = 0.25,
          min.pct = 0.1,
          only.pos = FALSE,
          test.use = "MAST",
          latent.vars = c("age", "sex"),
          verbose = FALSE
        )
        
        if (nrow(markers) > 0) {
          markers %>%
            rownames_to_column("gene") %>%
            mutate(
              sample = target_origin,
              celltype = ct,
              comparison_name = comp$name,
              case_group = paste(comp$case, collapse = "+"),
              control_group = paste(comp$control, collapse = "+"),
              case_n_cells = sum(seurat_ct$group_compare == "Case"),
              control_n_cells = sum(seurat_ct$group_compare == "Control"),
              case_n_patients = n_distinct(seurat_ct$orig.ident[seurat_ct$group_compare == "Case"]),
              control_n_patients = n_distinct(seurat_ct$orig.ident[seurat_ct$group_compare == "Control"])
            )
        } else {
          NULL
        }
      }, error = function(e) {
        message(sprintf("    [Error] cell type %s fail: %s", ct, conditionMessage(e)))
        return(NULL)
      })
    }) %>% compact()
    
    if (length(deg_results_list) > 0) {
      df_all_cells <- bind_rows(deg_results_list)
      
      df_filtered <- df_all_cells %>%
        filter(!str_detect(gene, "^RPS|^RPL|^MT-|^HLA-"))
      
      output_filename <- sprintf("Figure1_DEGs_%s_%s_fine.csv", target_origin, comp$name)
      write.csv(df_filtered, output_filename, row.names = FALSE)
  }
}

}

run_deg_pipeline(seurat_full = seurat_obj, target_origin = "Blood")
run_deg_pipeline(seurat_full = seurat_obj, target_origin = "CSF")

