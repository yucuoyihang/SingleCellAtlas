"""scTour trajectory inference for monocyte-state dynamics.

Negative-binomial modeling of UMI counts estimates monocyte pseudotime,
latent coordinates, and transcriptomic vector fields. A latent-space
neighborhood graph defines the UMAP representation.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import sctour as sct


N_TOP_GENES = 3000  # Maximum number of count-based Seurat v3 HVGs.
N_NEIGHBORS = 15  # Neighborhood size for the latent-space UMAP graph.
UMAP_MIN_DIST = 0.1  # Minimum distance controlling UMAP compactness.
RANDOM_SEED = 20260603  # Seed for model training and UMAP.


def use_raw_counts_in_x(adata: sc.AnnData) -> sc.AnnData:
    """Validate UMI counts and copy them to float32 .X for the scTour encoder."""
    if "counts" not in adata.layers:
        raise ValueError("scTour loss_mode='nb' requires raw UMI counts in adata.layers['counts'].")
    counts = adata.layers["counts"]
    values = counts.data if sp.issparse(counts) else np.asarray(counts)
    if (
        values.dtype.kind not in "iuf"
        or not np.all(np.isfinite(values))
        or np.any(values < 0)
        or np.any(values % 1 != 0)
    ):
        raise ValueError("adata.layers['counts'] must contain finite, nonnegative integers.")
    float_values = values.astype(np.float32)
    if not np.array_equal(values, float_values):
        raise ValueError("UMI counts must be exactly representable as float32 for scTour.")
    adata.X = (
        counts.tocsr().astype(np.float32, copy=True)
        if sp.issparse(counts)
        else np.array(counts, dtype=np.float32, copy=True)
    )
    return adata


def select_hvg_for_sctour(adata: sc.AnnData, n_top_genes: int) -> sc.AnnData:
    """Annotate QC before selecting up to n_top_genes count-based Seurat v3 HVGs.

    Retain pre-HVG n_genes_by_counts, which scTour uses for time orientation.
    This step annotates quality metrics; it does not filter cells.
    """
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    if adata.n_vars > n_top_genes:
        sc.pp.highly_variable_genes(
            adata,
            flavor="seurat_v3",
            n_top_genes=n_top_genes,
            subset=False,
        )
        hvg_mask = adata.var["highly_variable"].to_numpy(dtype=bool)
        adata = adata[:, hvg_mask].copy()
    else:
        adata.var["highly_variable"] = True
    return adata


def train_sctour_model(adata: sc.AnnData, train_percent: float | None):
    """Fit the NB model with equal encoder/ODE reconstruction weights (0.5).

    train_percent is the sampled cell fraction (0 < fraction <= 1); None
    retains scTour 1.0.0's default: 0.2 above 10,000 cells, 0.9 otherwise.
    Unspecified Trainer defaults: latent dimension 5; ODE/VAE hidden widths
    25/128; batch_norm=False; Euler solver; step_size=None (automatic choice);
    KL weight 1; batch size 1024; drop_last=False; Adam lr=0.001,
    wt_decay=1e-6, eps=0.01; validation fraction 0.1; use_gpu=True if available.
    nepoch=None gives min(round(400 * 10000 / round(n_cells * percent)), 400).
    random_state defaults to 0 upstream; here it is RANDOM_SEED.
    Source: https://sctour.readthedocs.io/en/stable/_modules/sctour/train.html
    """
    if train_percent is not None and (
        not np.isfinite(train_percent) or not 0 < train_percent <= 1
    ):
        raise ValueError("train_percent must be finite and in (0, 1], or None.")
    trainer_kwargs = dict(
        loss_mode="nb",
        alpha_recon_lec=0.5,
        alpha_recon_lode=0.5,
        random_state=RANDOM_SEED,
    )
    if train_percent is not None:
        trainer_kwargs["percent"] = train_percent
    tnode = sct.train.Trainer(adata, **trainer_kwargs)
    tnode.train()
    return tnode


def infer_sctour_dynamics(adata: sc.AnnData, tnode) -> sc.AnnData:
    """Store time, encoder/ODE latent states, their equal mixture, and derivatives.

    alpha_z=0.5 and alpha_predz=0.5 weight encoder and ODE latent coordinates,
    respectively. Evaluate the vector field at each cell's time and mixed state.
    """
    adata.obs["ptime"] = tnode.get_time()
    mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.5, alpha_predz=0.5)
    adata.obsm["X_TNODE"] = mix_zs
    adata.obsm["X_TNODE_z"] = zs
    adata.obsm["X_TNODE_pred_z"] = pred_zs
    adata.obsm["X_VF"] = tnode.get_vector_field(adata.obs["ptime"].to_numpy(), adata.obsm["X_TNODE"])
    return adata


def validate_root_fine(adata: sc.AnnData, root_fine: tuple[str, ...]) -> None:
    """Require every requested root label to occur among nonmissing fine labels."""
    if isinstance(root_fine, str):
        raise ValueError("root_fine must be a sequence of labels, not a single string.")
    if not root_fine:
        return
    if "fine" not in adata.obs:
        raise ValueError("root_fine requires adata.obs['fine'].")
    available = set(adata.obs["fine"].dropna().astype(str))
    missing = [label for label in root_fine if label not in available]
    if missing:
        raise ValueError(f"Unknown root_fine labels: {missing}")


def orient_pseudotime(adata: sc.AnnData, root_fine: tuple[str, ...]) -> sc.AnnData:
    """Apply an optional label-guided orientation heuristic, not biological evidence.

    Reverse time as max(ptime) - ptime and negate vectors only when the pooled
    root median exceeds the global median. This does not establish ancestry,
    guarantee that root cells are earliest, or refit the trained model.
    """
    validate_root_fine(adata, root_fine)
    if not root_fine:
        return adata
    root_mask = adata.obs["fine"].notna() & adata.obs["fine"].astype(str).isin(root_fine)
    root_median = float(np.nanmedian(adata.obs.loc[root_mask, "ptime"]))
    global_median = float(np.nanmedian(adata.obs["ptime"]))
    if np.isfinite(root_median) and root_median > global_median:
        adata.obs["ptime"] = np.nanmax(adata.obs["ptime"]) - adata.obs["ptime"]
        if "X_VF" in adata.obsm:
            adata.obsm["X_VF"] = -adata.obsm["X_VF"]
    return adata


def compute_umap_from_tnode(adata: sc.AnnData) -> sc.AnnData:
    """Build a neighborhood graph and UMAP embedding from the scTour latent space."""
    adata = adata[np.argsort(adata.obs["ptime"].to_numpy()), :].copy()
    sc.pp.neighbors(adata, use_rep="X_TNODE", n_neighbors=N_NEIGHBORS)
    sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, random_state=RANDOM_SEED)
    return adata


def plot_representative_sctour_panels(adata: sc.AnnData) -> dict[str, plt.Figure]:
    """Return one pseudotime UMAP and one vector-field figure without saving.

    adata is the annotated analysis result. Vector colors use fine labels when
    present, otherwise color=None. zs_key/vf_key identify latent coordinates
    and their derivatives; use_rep_neigh and t_key select latent/time neighbors.
    scTour's vector projection uses its default 20 neighbors and streamplot;
    size=80 and alpha=0.22 control cell marker area and opacity.
    The plotting API adds projection annotations to adata. Figures stay open.
    """
    if "ptime" not in adata.obs:
        raise ValueError("Plotting requires adata.obs['ptime'].")
    missing = [key for key in ("X_umap", "X_TNODE", "X_VF") if key not in adata.obsm]
    if missing:
        raise ValueError(f"Plotting requires adata.obsm keys: {missing}")
    color_key = "fine" if "fine" in adata.obs else None
    ptime_fig, ptime_ax = plt.subplots(figsize=(6.2, 5.4))
    sc.pl.umap(
        adata, color="ptime", frameon=False, show=False, save=False,
        cmap="viridis", ax=ptime_ax,
    )
    vector_fig, vector_ax = plt.subplots(figsize=(6.2, 5.4))
    sct.vf.plot_vector_field(
        adata,
        zs_key="X_TNODE",
        vf_key="X_VF",
        use_rep_neigh="X_TNODE",
        color=color_key,
        t_key="ptime",
        show=False,
        save=False,
        ax=vector_ax,
        legend_loc="none",
        frameon=False,
        size=80,
        alpha=0.22,
    )
    return {"pseudotime": ptime_fig, "vector_field": vector_fig}


def run_sctour_analysis(
    adata: sc.AnnData,
    root_fine: tuple[str, ...] = (),
    train_percent: float | None = None,
) -> tuple[sc.AnnData, sct.train.Trainer]:
    """Analyze an in-memory, final monocyte AnnData with layers['counts'].

    root_fine optionally identifies observed fine labels for heuristic time
    orientation; an empty tuple retains scTour's inferred orientation.
    train_percent overrides the training-cell fraction; None uses the Trainer
    defaults documented above. No normalization or cell filtering is applied.

    Return an annotated HVG-subset copy sorted by pseudotime, and the trained
    Trainer object. The input is unchanged. Trainer.adata retains training
    cell order, so align subsequent model outputs by obs_names, not position.
    Model time/vector predictions do not include the optional root heuristic.
    Plotting is separate and is never invoked automatically.
    """
    validate_root_fine(adata, root_fine)
    adata = use_raw_counts_in_x(adata.copy())
    adata = select_hvg_for_sctour(adata, N_TOP_GENES)
    tnode = train_sctour_model(adata, train_percent)
    adata = infer_sctour_dynamics(adata, tnode)
    adata = orient_pseudotime(adata, root_fine)
    adata = compute_umap_from_tnode(adata)

    return adata, tnode
