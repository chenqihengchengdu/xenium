# SPLIT: Spatial Purification of Layered Intracellular Transcripts

An R pipeline for processing 10X Xenium spatial transcriptomics data
with RCTD deconvolution and SPLIT purification.

Based on: Bilous, Buszta et al. *"From Transcripts to Cells: Dissecting
Sensitivity, Signal Contamination, and Specificity in Xenium Spatial
Transcriptomics."* bioRxiv (2025). doi:10.1101/2025.04.23.649965

## Workflow

```
Xenium raw data                 snRNA-seq reference
       │                              │
       │  cell_feature_matrix         │  ST16_Purity.rds
       │  cells.csv.gz                │  Anno_02_new.csv
       │                              │
       └──────────┬───────────────────┘
                  │
          ┌───────▼────────┐
          │  RCTD (spacexr) │   Cell type deconvolution
          │  doublet_mode   │   → proportions per spatial cell
          └───────┬────────┘
                  │
          ┌───────▼────────┐
          │  SPLIT purify   │   Signal purification
          │  linear model   │   → purified counts + metadata
          └───────┬────────┘
                  │
          ┌───────▼────────┐
          │  Seurat v4      │   SCTransform → PCA → UMAP
          │  cell annotation│   → marker gene validation
          └────────────────┘
```

*Optional: Proseg probabilistic segmentation can be run before RCTD
for improved transcript-to-cell assignment.*

## Requirements

### R packages

| Package  | Version | Source |
|----------|---------|--------|
| R        | ≥ 4.3   | CRAN |
| Seurat   | ≥ 4.4   | CRAN |
| spacexr  | ≥ 2.2   | GitHub: dmcable/spacexr |
| SPLIT    | ≥ 0.3   | GitHub: bdsc-tds/SPLIT |
| dplyr    | ≥ 1.1   | CRAN |
| ggplot2  | ≥ 3.5   | CRAN |
| tibble   | ≥ 3.2   | CRAN |

### Installation

```r
# Install from GitHub
remotes::install_github("dmcable/spacexr")
remotes::install_github("bdsc-tds/SPLIT")

# CRAN packages
install.packages(c("Seurat", "dplyr", "ggplot2", "tibble"))
```

### Optional: Proseg (probabilistic segmentation)

```bash
# Via conda (recommended)
conda create -n proseg -c bioconda -c conda-forge rust-proseg

# Or via cargo
cargo install proseg
```

Proseg: Jones et al., *Nature Methods* (2025). https://github.com/dcjones/proseg

## Input Data

### Xenium spatial data (per sample)
```
sample_output/
├── cell_feature_matrix/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── cells.csv.gz              # x_centroid, y_centroid, cell_id
└── transcripts.csv.gz         # (optional, for Proseg)
```

### snRNA-seq reference
- `reference.rds` — Seurat object with RNA assay (counts)
- `annotations.csv` — cell type labels

## Usage

### 1. (Optional) Proseg segmentation
```bash
bash 01_proseg_segmentation.sh
```

### 2. RCTD deconvolution + SPLIT purification
```r
# Edit configuration at top of script:
#   REF_RDS, ANNO_CSV, DATA_BASE, SAMPLE_PATTERN
Rscript 02_split_purify.R
```

### 3. Cell type annotation
```r
Rscript 03_cell_annotation.R
```

## Output

```
project/
├── data/
│   ├── purified_xe_output-<sample_id>.rds    # purified Seurat objects
│   └── umap_<sample_id>.png                   # UMAP previews
├── output/
│   ├── <sample_id>_01_umap_first_type.png     # annotated UMAP
│   ├── <sample_id>_02_markers.png             # marker gene FeaturePlots
│   ├── <sample_id>_03_dotplot.png             # marker DotPlot
│   ├── <sample_id>_04_celltype_proportions.csv # cell type stats
│   └── <sample_id>_05_epi_subtypes.png        # epithelial sub-clusters
└── proseg_output/                              # (optional) Proseg results
```

## Important Notes

1. **max_cores = 1**: Use single-core RCTD to avoid parallel worker deadlock
   with large datasets. Set in `create.RCTD(puck, ref, max_cores = 1)`.

2. **Save before plotting**: `saveRDS()` should precede `ggsave()` so data
   is preserved even if visualization fails.

3. **Timing**: RCTD+SPLIT takes 5–10 hours per sample (70K–115K cells,
   18 cell types) with single-core execution.

4. **Memory**: Expect 30–50 GB RSS per sample during processing.

5. **Reference downsampling**: spacexr downsamples to 10,000 cells per type.
   This is expected behavior and does not affect annotation quality.

## Citation

If you use this pipeline, please cite:

- **SPLIT method**: Bilous M, Buszta D, Bac J et al. From Transcripts to
  Cells: Dissecting Sensitivity, Signal Contamination, and Specificity in
  Xenium Spatial Transcriptomics. *bioRxiv* (2025).

- **RCTD**: Cable DM, Murray E, Zou LS et al. Robust decomposition of cell
  type mixtures in spatial transcriptomics. *Nature Biotechnology*
  40:517–526 (2022).

- **Seurat**: Hao Y et al. Dictionary learning for integrative, multimodal
  and scalable single-cell analysis. *Nature Biotechnology* 42:293–304 (2024).

## License

This pipeline is adapted from the analysis scripts of the SPLIT paper
authors (huangl lab). Original code at:
https://github.com/bdsc-tds/xenium_analysis_pipeline
