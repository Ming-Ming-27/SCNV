from anndata import AnnData

def subset_shared_genes(
    adata: AnnData,
    gene_pos_df
):
    """
    Subset adata and gene_pos_df to their shared genes.
    """
    shared = adata.var_names.intersection(gene_pos_df.index)
    adata._inplace_subset_var(shared)
    gene_pos_df = gene_pos_df.loc[shared].copy()
    return adata, gene_pos_df