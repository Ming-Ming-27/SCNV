import numpy as np
import pandas as pd
from scipy import sparse
import warnings
from .subset import subset_shared_genes

def infer_cnv_and_annotate(
    adata,
    gene_pos_df,
    celltype_key="cell_type",
    normal_key="normal_reference",
    window_size=100,
    log1p=True,
    thresh_del=-1.0,
    thresh_half=0.0,
    thresh_gain=1.0,
    min_cells=10,
    min_annots=5,
    layer="counts"
):
    """
    Infers CNVs (copy number variations) based on sliding window expression differences,
    and writes the results to `adata.obs['annotated_cnvs']`.

    - For each `cell_type`, the baseline is computed using cells labeled as 'normal' in `normal_reference`
    - The window size for chromosomes X/Y is reduced to `window_size // 5` (minimum 1)
    - Only segments that appear in at least `min_cells` cells are kept
    - Only full annotation strings that occur in at least `min_annots` cells are retained
    """
    # 1. Extract raw expression matrix
    data = adata.layers.get(layer, adata.X)
    mat = data.toarray() if sparse.issparse(data) else data
    expr = pd.DataFrame(mat, index=adata.obs_names, columns=adata.var_names)

    # 2. Apply log1p transformation
    if log1p:
        expr = np.log1p(expr)

    # 3. Build sliding windows and compute mean expression
    windows, labels = [], []
    for chrom, grp in gene_pos_df.groupby("chromosome"):
        grp_sorted = grp.sort_values("start")
        genes = grp_sorted.index.tolist()
        # Use smaller window size for sex chromosomes
        ws = window_size // 5 if chrom in ('X') else window_size
        ws = max(1, ws)
        for i in range(0, len(genes), ws):
            w = genes[i:i+ws]
            if not w:
                continue
            windows.append(w)
            s = grp_sorted.loc[w, "start"].min()
            e = grp_sorted.loc[w, "end"].max()
            labels.append(f"{chrom}:{int(s)}-{int(e)}")

    # Calculate per-cell mean expression for each window
    win_df = pd.DataFrame(
        {lbl: expr[w].mean(axis=1) for lbl, w in zip(labels, windows)},
        index=expr.index
    )

    # 4. Compute baseline expression using normal cells in each cell type
    ct = adata.obs[celltype_key]
    is_ref = adata.obs[normal_key] == 'normal'
    baseline_df = win_df[is_ref].groupby(ct[is_ref]).mean()
    baseline_for_cells = baseline_df.loc[ct.values]
    baseline_for_cells.index = ct.index

    # 5. Compute deviation matrix
    cnv_df = win_df.subtract(baseline_for_cells, axis=0)

    # 6. Make raw CNV calls for each cell
    raw_annots = []
    for _, row in cnv_df.iterrows():
        segs = []
        for lbl, val in row.items():
            if val < thresh_del:
                call = 0
            elif val < thresh_half:
                call = 1
            elif val > thresh_gain:
                call = 4
            else:
                continue
            segs.append(f"{lbl} (CN {call})")
        raw_annots.append(", ".join(segs))

    # 7. Filter by segment frequency (keep segments seen in ≥ min_cells)
    all_segs = [seg.strip() for s in raw_annots if s for seg in s.split(",")]
    seg_counts = pd.Series(all_segs).value_counts()
    allowed = set(seg_counts[seg_counts >= min_cells].index)

    filtered = []
    for s in raw_annots:
        if not s:
            filtered.append("")
            continue
        kept = [seg.strip() for seg in s.split(",") if seg.strip() in allowed]
        filtered.append(", ".join(kept))

    adata.obs['annotated_cnvs'] = filtered

    # 8. Keep only annotation patterns that occur in ≥ min_annots cells
    vc = adata.obs['annotated_cnvs'].value_counts()
    keep = set(vc[vc >= min_annots].index)
    adata.obs['annotated_cnvs'] = adata.obs['annotated_cnvs'].where(
        adata.obs['annotated_cnvs'].isin(keep), ""
    )

    return cnv_df



def annotate_cnv(
    adata,
    gene_pos_df,
    celltype_key="cell_type",
    window_size=50,
    thresh_del=-0.15,
    thresh_half=-0.1,
    thresh_gain=0.3,
    min_cells=400,
    min_annots=50,
    layer="counts"
):
    """
    Runs CNV inference chromosome-by-chromosome and merges results into
    `adata.obs['annotated_cnvs']` as a comma-separated annotation string.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix, must contain .var['chromosome'] field
    gene_pos_df : pd.DataFrame
        Gene position DataFrame indexed by gene name, must contain 'chromosome', 'start', 'end'
    celltype_key : str
        Key for cell type annotation in `adata.obs`
    window_size : int
        Number of genes per sliding window
    thresh_del, thresh_half, thresh_gain : float
        Thresholds for calling CNV states
    min_cells : int
        Keep only segments found in at least this many cells
    min_annots : int
        Keep only annotations occurring at least this many times
    layer : str
        Expression layer to use (e.g., 'counts')

    Returns
    -------
    adata : AnnData
        Updated `adata` object with `obs['annotated_cnvs']` column populated
    """
    # Initialize (or reset) the annotation field
    adata.obs['annotated_cnvs'] = [''] * adata.n_obs

    # Iterate through chromosomes
    for chrom in adata.var['chromosome'].unique():
        # Subset adata and gene position table for this chromosome
        adchr = adata[:, adata.var['chromosome'] == chrom].copy()
        adchr, gene_pos_df_sub = subset_shared_genes(adchr, gene_pos_df)

        # Suppress warnings and run CNV inference
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=FutureWarning)
            infer_cnv_and_annotate(
                adchr,
                gene_pos_df_sub,
                celltype_key=celltype_key,
                window_size=window_size,
                thresh_del=thresh_del,
                thresh_half=thresh_half,
                thresh_gain=thresh_gain,
                min_cells=min_cells,
                min_annots=min_annots,
                layer=layer
            )

        # Merge annotations back into main adata
        for cell_id, annotation in adchr.obs['annotated_cnvs'].items():
            if not annotation:
                continue
            existing = adata.obs.at[cell_id, 'annotated_cnvs']
            if existing:
                adata.obs.at[cell_id, 'annotated_cnvs'] = f"{existing}, {annotation}"
            else:
                adata.obs.at[cell_id, 'annotated_cnvs'] = annotation

    return adata
