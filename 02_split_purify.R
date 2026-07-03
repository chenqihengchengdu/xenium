# ================================================================
# 02_split_purify.R
# SPLIT Spatial Purification Pipeline for Xenium Data
#
# Core workflow: RCTD deconvolution → SPLIT purification → Seurat
#
# Reference:
#   Bilous, Buszta et al. "From Transcripts to Cells."
#   bioRxiv (2025). doi:10.1101/2025.04.23.649965
#
# Dependencies:
#   R >= 4.3, Seurat >= 4.4, spacexr >= 2.2, SPLIT >= 0.3
#   GitHub: dmcable/spacexr, bdsc-tds/SPLIT
#
# Input:
#   - Xenium cell_feature_matrix + cells.csv.gz (per sample)
#   - snRNA-seq reference Seurat object (.rds)
#   - Reference cell type annotations (.csv)
#
# Output:
#   - data/purified_xe_output-<sample_id>.rds
#   - data/umap_<sample_id>.png
# ================================================================

library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)
library(spacexr)
library(SPLIT)

# ======================== Configuration ========================
# --- File paths ---
REF_RDS   <- "/path/to/ST16_Purity.rds"       # snRNA-seq reference
ANNO_CSV  <- "/path/to/Anno_02_new.csv"       # reference cell type labels
DATA_BASE <- "/path/to/xenium/data"           # Xenium output directory
OUTPUT_DIR <- "data"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Sample filter (regex) ---
# "P020|P021" for testing; "P0" for all samples
SAMPLE_PATTERN <- "P020|P021"

# ======================== Load Reference ========================
message("Loading snRNA-seq reference...")
ref_srt <- readRDS(REF_RDS)
anno    <- read.csv(ANNO_CSV, row.names = 1)
ref_srt <- AddMetaData(ref_srt, anno)

counts_ref     <- ref_srt@assays$RNA$counts
cell_types_ref <- ref_srt$Anno_02_new %>% as.factor()
reference      <- Reference(counts = counts_ref, cell_types = cell_types_ref)

message(sprintf("Reference: %d cells, %d cell types",
                ncol(counts_ref), length(levels(cell_types_ref))))

# ======================== Core Function ========================
# Adapted from split_fun.r (huangl lab)
#
# Processes one Xenium sample through the full pipeline:
#   Read10X → RCTD → SPLIT post-process → SPLIT purify → Seurat
#
# @param raw_data_path Path to Xenium cell_feature_matrix directory
# @param ref           spacexr Reference object (snRNA-seq)
# @param id            Sample identifier string
# @param output_dir    Directory for output RDS and PNG files
# @return              Purified Seurat object
process_raw_data_with_rctd <- function(raw_data_path, ref, id, output_dir) {

    # ---------- 1. Read Xenium data ----------
    raw    <- Read10X(raw_data_path)
    counts <- raw[[1]]
    message(sprintf("[%s] %d genes × %d cells", id, nrow(counts), ncol(counts)))

    # ---------- 2. Read spatial coordinates ----------
    coords_path <- sub("cell_feature_matrix$", "cells.csv.gz", raw_data_path)
    coords <- read.csv(coords_path) %>%
        select(cell_id, x_centroid, y_centroid) %>%
        column_to_rownames(var = "cell_id") %>%
        rename(xcoord = x_centroid, ycoord = y_centroid)

    # ---------- 3. Build SpatialRNA puck ----------
    puck <- SpatialRNA(coords = coords, counts = counts)

    # ---------- 4. RCTD deconvolution ----------
    # max_cores = 1 is CRITICAL to avoid parallel worker deadlock
    message(sprintf("[%s] Running RCTD (doublet mode)...", id))
    myRCTD <- create.RCTD(puck, ref, max_cores = 1)
    myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")

    # ---------- 5. SPLIT post-processing ----------
    message(sprintf("[%s] SPLIT post-process RCTD...", id))
    myRCTD <- SPLIT::run_post_process_RCTD(myRCTD)

    # ---------- 6. SPLIT purification ----------
    # Purifies the count matrix by decomposing mixed transcript signals
    # using a linear model informed by RCTD cell type proportions.
    # DO_purify_singlets = TRUE also refines single-cell signals.
    message(sprintf("[%s] SPLIT purify...", id))
    res_split <- SPLIT::purify(
        counts             = counts,
        rctd               = myRCTD,
        DO_purify_singlets = TRUE
    )

    # ---------- 7. Build purified Seurat object ----------
    xe_purified <- CreateSeuratObject(
        counts    = res_split$purified_counts,
        meta.data = res_split$cell_meta,
        assay     = "Xenium"
    )

    n_before <- ncol(counts)
    n_after  <- ncol(xe_purified)
    message(sprintf("[%s] Purification: %d → %d cells", id, n_before, n_after))

    # ---------- 8. Filter, normalize, reduce dimensions ----------
    xe_purified <- subset(xe_purified, subset = nCount_Xenium > 5) %>%
        SCTransform(assay = "Xenium") %>%
        RunPCA() %>%
        RunUMAP(dims = 1:30)

    # ---------- 9. Save purified data (before plotting!) ----------
    rds_path <- file.path(output_dir, paste0("purified_xe_output-", id, ".rds"))
    saveRDS(xe_purified, file = rds_path)
    message(sprintf("[%s] Saved: %s", id, rds_path))

    # ---------- 10. UMAP visualization ----------
    p <- UMAPPlot(xe_purified, group.by = c("first_type"),
                  label = TRUE, repel = TRUE) + theme(aspect.ratio = 1)
    ggsave(file.path(output_dir, paste0("umap_", id, ".png")),
           plot = p, width = 8, height = 6, dpi = 150)

    return(xe_purified)
}

# ======================== Discover & Filter Samples ========================
all_samples <- list.dirs(DATA_BASE, full.names = TRUE, recursive = TRUE) %>%
    grep(pattern = "cell_feature_matrix$", value = TRUE)

selected <- grep(SAMPLE_PATTERN, all_samples, value = TRUE)
message(sprintf("Found %d samples, selected %d for processing",
                length(all_samples), length(selected)))

# ======================== Batch Processing ========================
results <- list()
for (i in seq_along(selected)) {
    raw_path  <- selected[i]
    sample_id <- basename(dirname(raw_path))

    message(sprintf("\n======== [%d/%d] %s ========",
                    i, length(selected), sample_id))
    t_start <- Sys.time()

    tryCatch({
        res <- process_raw_data_with_rctd(
            raw_data_path = raw_path,
            ref           = reference,
            id            = sample_id,
            output_dir    = OUTPUT_DIR
        )
        results[[sample_id]] <- res
        t_elapsed <- round(difftime(Sys.time(), t_start, units = "hours"), 1)
        message(sprintf("[%s] Complete (%.1f h)", sample_id, t_elapsed))
    }, error = function(e) {
        message(sprintf("[%s] FAILED: %s", sample_id, e$message))
    })
}

message(sprintf("\n=== Done: %d/%d samples succeeded ===",
                length(results), length(selected)))
