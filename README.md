This repository contains analysis scripts for investigating single-cell transcriptomic atlas of Alzheimer’s disease (AD). The analyses include single-cell RNA-sequencing of blood and cerebrospinal fluid (CSF), monocyte-state characterization, cell–cell communication, trajectory inference, validation using external data, and genetic and proteomic profiling.

The scripts are organized as research analysis modules and figure-generation workflows. They are not packaged as a single, fully automated pipeline; most scripts expect study-specific input files, metadata, and analysis objects.

## Repository contents

| Script | Description |
| --- | --- |
| `01_UMAP_and_Marker_gene.R` | Visualize annotated cell populations and marker genes. |
| `02_cell_proportion_calculation.R` | Calculate major cell-type and fine-subtype proportions. |
| `03_cell_proportion_significance.R` | Test changes in cell proportions across AD diagnostic groups. |
| `04_DEG_calculation.R` | Perform differential gene-expression analyses across diagnostic comparisons. |
| `05_DEG_cluster.R` | Cluster and analyze differential-expression patterns across adjacent AD stages. |
| `06_Monocyte_analysis.R` | Analyze monocyte populations and generate associated figures and enrichment results. |
| `07_brain_plot.m` | Visualize cortical and subcortical results using the ENIGMA MATLAB toolbox. |
| `08_External_data_validation.R` | Validate selected cell and gene-expression findings in external datasets. |
| `09_Monocle2_trajectory_BEAM_and_visualization.R` | Infer Monocle 2 trajectories and analyze branch-dependent expression. |
| `10_sctour_monocyte_trajectory.py` | Infer monocyte trajectories and latent dynamics using scTour. |
| `11_Cellfate_signature_scoring_and_visualization.R` | Score and visualize cell-fate gene signatures. |
| `12_Cellchat_analysis_and_visualization.R` | Analyze cell–cell communication, including SPP1 signaling. |
| `13_Peripheral_and central_analysis.R` | Explore the relationships of peripheral and central immune cell features.
| `14_scDRS.sh` | Run scDRS scoring and downstream cell-type association analyses. |
| `15_Differential_protein_analysis.R` | Test protein associations with diagnostic groups. |
| `16_Protein_abundance_external_validation.R` | Visualize protein abundance and validate findings in external cohorts. |

## Requirements

The analyses use R, Python, and MATLAB. Required software and packages vary by script.

- **R:** Packages used across the R scripts include Seurat, tidyverse, MAST, emmeans, Monocle 2, AUCell, GSEABase, CellChat, edgeR, and data.table. Some analyses also use additional plotting, enrichment, or statistical packages.
- **Python:** The scTour workflow uses Python packages including Scanpy, scTour, NumPy, SciPy, and Matplotlib.
- **Command line:** The scDRS workflow requires scDRS and a Conda/Mamba environment.
- **MATLAB:** The brain visualization script requires the ENIGMA toolbox and its MATLAB dependencies.

Install the software and package versions required by the specific analyses before running them. The repository does not currently provide a single environment lockfile.

## Data and inputs

The scripts refer to study-specific files that are not included in this repository, such as Seurat `.rds` objects, count matrices, metadata tables, differential-expression results, protein-abundance data, and external validation datasets. Some analyses also require reference resources, including ligand–receptor networks and gene-set files.

Obtain the relevant data and resources from their original sources, following their access conditions and licenses. Update file paths and, where necessary, column names to match your local data before running a script. Do not upload restricted or identifiable participant-level data.

## Running the analyses

1. Obtain the input data and reference resources required by the analysis of interest.
2. Install the corresponding language runtime, packages, and external tools.
3. Review the input paths, expected metadata fields, and analysis parameters in the script.
4. Run the script in the appropriate R, Python, shell, or MATLAB environment.

Scripts are intended to be run independently or as parts of a larger analysis workflow; the numeric prefixes indicate an approximate analysis order, not a guaranteed end-to-end execution sequence. Outputs are generally written to the working directory or to paths defined in each script.

## Reproducibility

Several scripts set random seeds for stochastic analyses. Results may still vary with software versions, input data, hardware, and package defaults. For reproducibility, record the versions of the software and packages used, along with any local changes to paths or parameters.

## Citation

If you use these scripts, please cite the associated study and the methods or software packages used. 
