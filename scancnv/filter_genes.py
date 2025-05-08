import re
from anndata import AnnData

def filter_genes_by_chromosome(
    adata: AnnData,
    chrom_key: str = 'chromosome',
    inplace: bool = False
) -> AnnData:
    """
    Keep only genes on chromosomes 1-22, X, Y.
    """
    chroms = adata.var[chrom_key].astype(str)
    pattern = r'^(?:[1-9]|1[0-9]|2[0-2]|X|Y)$'
    mask = chroms.str.match(pattern)

    if inplace:
        adata._inplace_subset_var(mask.values)
        return adata
    else:
        return adata[:, mask.values].copy()