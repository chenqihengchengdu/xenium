#!/bin/bash
# ============================================================
# 01_proseg_segmentation.sh
# Run Proseg probabilistic cell segmentation on Xenium data
#
# Proseg v3.1.1 (Jones et al., Nature Methods, 2025)
# GitHub: https://github.com/dcjones/proseg
#
# Input:  Xenium transcripts.csv.gz (per sample)
# Output: proseg_output/<sample_id>/
#         ├── expected-counts.csv.gz    # non-integer expected counts
#         ├── cell-metadata.csv.gz      # cell centroids, volume, etc.
#         ├── counts.csv.gz             # integer point-estimate counts
#         └── proseg-output.zarr/       # full spatialdata output
#
# NOTE: Proseg is OPTIONAL. The main SPLIT pipeline works directly
#       with Xenium cell_feature_matrix output.
# ============================================================

set -e

PROSEG="/path/to/proseg"
DATA_BASE="/path/to/xenium/data"
OUTPUT_BASE="./proseg_output"
NTHREADS=20

echo "=== Proseg: Probabilistic Segmentation of Xenium Data ==="
echo "Data:   $DATA_BASE"
echo "Output: $OUTPUT_BASE"
echo "Threads: $NTHREADS"
echo ""

# Iterate over all Xenium sample directories
for sample_dir in "$DATA_BASE"/output-XETG00365__*; do
    sample_id=$(basename "$sample_dir")
    transcript_csv="$sample_dir/transcripts.csv.gz"

    if [ ! -f "$transcript_csv" ]; then
        echo "[SKIP] $sample_id — no transcripts.csv.gz"
        continue
    fi

    sample_out="$OUTPUT_BASE/$sample_id"
    mkdir -p "$sample_out"

    echo "========================================"
    echo "[RUN] $sample_id"
    echo "  Input:  $transcript_csv"
    echo "  Output: $sample_out"
    echo ""

    # Proseg with Xenium preset
    #   --samples 200:           MCMC main sampling iterations
    #   --burnin-samples 200:    coarse burn-in iterations
    #   --nthreads N:            CPU threads for parallelization
    #   --output-expected-counts:  non-integer expected expression matrix
    #   --output-cell-metadata:    cell-level metadata
    #   --output-counts:           integer count matrix
    #   --output-spatialdata:      full spatialdata zarr archive
    "$PROSEG" \
        --xenium \
        --samples 200 \
        --burnin-samples 200 \
        --nthreads "$NTHREADS" \
        --output-expected-counts "$sample_out/expected-counts.csv.gz" \
        --output-cell-metadata "$sample_out/cell-metadata.csv.gz" \
        --output-counts "$sample_out/counts.csv.gz" \
        --output-spatialdata "$sample_out/proseg-output.zarr" \
        "$transcript_csv"

    echo "[DONE] $sample_id"
    echo ""
done

echo "=== All samples complete ==="
