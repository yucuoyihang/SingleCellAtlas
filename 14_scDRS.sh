# scDRS analysis ===============================================================

# Activate the scDRS environment
mamba activate scdrs_env

# Input and output files
H5AD_FILE="ADCN_allfine_group.h5ad"
GS_FILE="processed_geneset_human.gs"
OUT_DIR="scdrs_results"

# 1. Compute cell-level scDRS scores
scdrs compute-score \
  --h5ad-file "${H5AD_FILE}" \
  --h5ad-species human \
  --gs-file "${GS_FILE}" \
  --gs-species human \
  --out-folder "${OUT_DIR}" \
  --cov-file None \
  --flag-filter-data True \
  --flag-raw-count True \
  --n-ctrl 1000 \
  --flag-return-ctrl-raw-score False \
  --flag-return-ctrl-norm-score True

# 2. Test trait associations across CSF and blood fine celltypes
scdrs perform-downstream \
  --h5ad-file "${H5AD_FILE}" \
  --score-file "${OUT_DIR}/@.full_score.gz" \
  --out-folder "${OUT_DIR}" \
  --group-analysis origin_finecelltype \
  --flag-filter-data True \
  --flag-raw-count True