# scancnv

**scancnv** is a lightweight toolkit for **copy‑number variation (CNV) inference, annotation, and evaluation** from single‑cell RNA‑seq count data.  
It implements a sliding‑window expression–difference strategy with automatic diploid‑reference selection, making it suitable for datasets that lack genome‑wide coverage or external CNV references.

| Module / Function            | Purpose                                                         |
| ---------------------------- | --------------------------------------------------------------- |
| `filter_genes_by_chromosome` | Keep standard chromosomes 1‑22 + X/Y                            |
| `subset_shared_genes`        | Synchronise genes between an `AnnData` object and position table|
| `annotate_normal_cells`      | Mark diploid reference cells (per cell type)                    |
| `infer_cnv_and_annotate`     | Sliding‑window CNV calling & per‑cell annotation                |
| `annotate_cnv`               | Chromosome‑wise CNV calling and global merge                    |
| `assess_cnv`                 | Compare predictions with ground‑truth labels                    |

---

## Installation

```bash
# 1 – create a fresh environment (recommended)
conda create -n scancnv python=3.10 anndata pandas numpy scipy scanpy
conda activate scancnv

# 2 – install from source
git clone https://github.com/YourLab/scancnv.git
cd scancnv
pip install -e .
```

### Dependencies

| Package                    | Version ≥ | Role                     |
| -------------------------- | --------- | ------------------------ |
| `anndata`                  | 0.9       | single‑cell container    |
| `pandas`, `numpy`, `scipy` | —         | data wrangling           |
| `scanpy`                   | optional  | I/O & plotting helpers   |

---

## Quick Start

```python
import scancnv as cnv
import scanpy as sc

# 1. Load AnnData
adata = sc.read_h5ad("pbmc_10k_v3_filtered_feature_bc_matrix.h5ad")

# 2. Keep standard chromosomes only
adata = cnv.filter_genes_by_chromosome(adata, inplace=True)

# 3. Prepare a gene‑position table
gene_pos = adata.var[["chromosome", "start", "end"]]

# 4. Mark diploid reference cells
cnv.annotate_normal_cells(
    adata,
    mode="inferred",       # or "simulated"
    min_refs=100
)

# 5. Chromosome‑wise CNV inference
cnv.annotate_cnv(
    adata,
    gene_pos_df=gene_pos,
    window_size=100,
    thresh_del=-0.8,
    thresh_half=-0.2,
    thresh_gain=0.5,
    min_cells=20,
    min_annots=5,
    layer="counts"
)

# 6. (Optional) Evaluate against simulated truth
cnv.assess_cnv(adata)
```

---
## DATASET
please download files from googledrive: https://drive.google.com/drive/folders/1gXJyu3XWCh6ZjN99JNyX-N5gW9ts4jkQ?usp=drive_link.
## Data Requirements

| Object                           | Mandatory fields                                    | Notes                                                         |
| -------------------------------- | --------------------------------------------------- | ------------------------------------------------------------- |
| `AnnData.X` / `layers['counts']` | raw count matrix                                   | genes × cells                                                 |
| `adata.var`                      | `chromosome`, `start`, `end`                       | genomic coordinates                                           |
| `adata.obs`                      | `cell_type`*                                       | cell‑type labels (string)                                     |
| `gene_pos_df`                    | index = gene name; cols `chromosome`,`start`,`end` | if already present in `adata.var`, simply pass `adata.var`    |

\* Column names are configurable via function parameters.

---

## Core API

### `filter_genes_by_chromosome`

```python
filter_genes_by_chromosome(adata, chrom_key="chromosome", inplace=False)
```

Removes non‑standard/mitochondrial genes, retaining chr 1‑22 + X/Y.

---

### `annotate_normal_cells`

```python
annotate_normal_cells(
    adata,
    simulated_key="simulated_cnvs",
    inferred_key="annotated_cnvs",
    celltype_key="cell_type",
    min_refs=50,
    mode="simulated" | "inferred",
    col_name="normal_reference"
)
```

| Mode        | Logic                                                                           |
| ----------- | ------------------------------------------------------------------------------- |
| `simulated` | Cells whose `simulated_cnvs` field is **empty** are marked as *normal*          |
| `inferred`  | For each `cell_type`, pick `min_refs` cells with the **fewest inferred segments* |

---

### `infer_cnv_and_annotate`

```python
infer_cnv_and_annotate(
    adata_sub,
    gene_pos_df_sub,
    window_size=100,
    log1p=True,
    thresh_del=-1.0,
    thresh_half=0.0,
    thresh_gain=1.0,
    min_cells=10,
    min_annots=5,
    layer="counts"
)
```

| Step       | Description                                                                              |
| ---------- | ---------------------------------------------------------------------------------------- |
| Windowing  | Fixed **N genes** per window; window size for X/Y is automatically divided by 5          |
| Baseline   | Mean expression of diploid references within the same `cell_type`                        |
| Calling    | Three‑threshold rule: < `thresh_del`→CN0; `thresh_del`–`thresh_half`→CN1; > `thresh_gain`→CN4 |
| Filtering  | Keep segments present in ≥ `min_cells` cells and annotation patterns appearing ≥ `min_annots` times |

---

### `annotate_cnv`

Runs `infer_cnv_and_annotate` **per chromosome** and concatenates the resulting strings into `obs['annotated_cnvs']`.

---

### `assess_cnv`

```python
assess_cnv(
    adata,
    pred_key="annotated_cnvs",
    label_key="simulated_cnvs"
)
```

Per‑cell × per‑chromosome comparison.  
A chromosome counts as **true positive** if any predicted segment overlaps a true segment **and** the copy‑number is identical.

Outputs TP, FP, FN, TN and overall **Accuracy**, **Precision**, **Sensitivity**.

---

## Parameter‑Tuning Tips

| Scenario                   | Suggestion                                       |
| -------------------------- | ----------------------------------------------- |
| High noise / low depth     | Decrease `window_size` and/or relax thresholds  |
| Few diploid cells          | Increase `min_refs`; ensure a robust baseline   |
| TPM / CPM matrices         | Disable `log1p` or adjust thresholds accordingly|

---

## FAQ

| Question                     | Resolution                                                     |
| ---------------------------- | -------------------------------------------------------------- |
| `KeyError: 'chromosome'`     | Ensure `adata.var` contains a `chromosome` column              |
| Few or zero CNV calls        | Lower `min_cells` / `min_annots` or inspect threshold settings|
| Evaluation metrics are `nan` | Provide ground‑truth labels in `label_key`                    |

---

## Citation

---

## License

Released under the **MIT License**.

---

## Maintainer
