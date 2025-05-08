import pandas as pd
from anndata import AnnData

def annotate_normal_cells(
    adata: AnnData,
    simulated_key: str = 'simulated_cnvs',
    inferred_key: str = 'annotated_cnvs',
    celltype_key: str = 'cell_type',
    min_refs: int = 50,
    mode: str = 'simulated',
    col_name: str = 'normal_reference'
):
    """
    Add a new column `col_name` to `adata.obs`:
      - Cells selected as diploid normal references will be labeled as 'normal'
      - All other cells will be labeled as an empty string

    The selection strategy is determined by the `mode` parameter:
      * 'simulated': all cells with an empty `simulated_key` value are marked as normal (no fallback logic)
      * 'inferred': for each cell type, select the top `min_refs` cells with the fewest CNV segments (based on `inferred_key`)

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    simulated_key : str
        Column name in `obs` that contains simulated ground truth CNV labels
    inferred_key : str
        Column name in `obs` containing inferred CNV labels
    celltype_key : str
        Column name indicating the cell type
    min_refs : int
        Number of normal reference cells to select per cell type in 'inferred' mode
    mode : {'simulated', 'inferred'}
        Strategy for selecting normal reference cells
    col_name : str
        Name of the new column to be added in `adata.obs`
    """
    obs = adata.obs.copy()

    # Helper function to count the number of CNV segments in a label string
    def count_segs(s):
        if pd.isna(s) or not s.strip():
            return 0
        return len([seg for seg in s.split(',') if seg.strip()])

    # Initialize the new column with empty strings
    obs[col_name] = ''

    # Process each cell type separately
    for ct, sub in obs.groupby(celltype_key):
        if mode == 'simulated':
            # In simulated mode, mark cells with empty CNV labels as normal
            normals = sub.index[sub[simulated_key].fillna('') == ''].tolist()
            chosen = normals

        else:  # mode == 'inferred'
            # In inferred mode, select the top `min_refs` cells with the fewest CNV segments
            seg_counts = sub[inferred_key].fillna('').apply(count_segs)
            sorted_cells = seg_counts.sort_values().index.tolist()
            chosen = sorted_cells[:min_refs]

        # Mark selected cells as normal
        obs.loc[chosen, col_name] = 'normal'

    # Write the result back to `adata.obs`
    adata.obs[col_name] = obs[col_name]
