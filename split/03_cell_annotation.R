# ================================================================
# 03_cell_annotation.R
# Cell type annotation and visualization for SPLIT-purified data
#
# Reads purified_xe_output-*.rds from the data/ directory and
# generates UMAP plots, FeaturePlots, DotPlots, and cell type
# proportion tables.
#
# Marker genes are curated from the PCa (prostate cancer) literature
# and the Xenium SPLIT paper (Bilous, Buszta et al., 2025).
# ================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# ======================== Configuration ========================
DATA_DIR   <- "data"
OUTPUT_DIR <- "output"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ======================== Marker Gene Panels ========================
# Hierarchical cell type markers covering major lineages detected
# by Xenium spatial transcriptomics with SPLIT purification.

MARKER_PANELS <- list(

    # --- Epithelial lineage ---
    "Epithelial_Luminal" = c(
        "TMPRSS2", "NKX3-1", "KLK3", "PCA3",    # prostate luminal
        "ACPP", "CDH1", "AR", "EPCAM"            # pan-epithelial
    ),
    "Epithelial_Club" = c(
        "OLFM4", "KRT13", "SOX9"                 # club / hillock
    ),
    "Epithelial_Basal" = c(
        "FGFR2", "TP63", "KRT5", "CSTB", "LYPD3" # basal / stem
    ),

    # --- Stromal lineage ---
    "Fibroblast" = c(
        "COL1A1", "PDPN", "FAP", "PDGFRB",
        "IGF1", "CCDC80", "SFRP1", "DCN", "LUM"
    ),
    "SMC_Pericyte" = c(
        "ACTA2", "MYLK", "MYH11", "TAGLN",
        "RGS5", "CSPG4", "PDGFRB"
    ),

    # --- Endothelial lineage ---
    "Endothelial" = c(
        "PECAM1", "VWF", "ACKR1", "PLVAP",
        "A2M", "CD34", "CLDN5", "ENG"
    ),

    # --- Immune: Myeloid lineage ---
    "Myeloid" = c(
        "CD68", "CD163", "CD14", "CSF1R",
        "HLA-DRA", "HLA-DPA1", "HLA-DPB1",
        "CXCL8", "C1QB", "CD1C"
    ),

    # --- Immune: Lymphoid lineage ---
    "T_Cell" = c(
        "CD3E", "CD8A", "CD4", "CCL5",
        "CD160", "NKG7", "GNLY", "CD247",
        "IL7R", "FOXP3"
    ),
    "B_Cell" = c(
        "MS4A1", "CD79A", "CD79B", "CXCR5",
        "BANK1", "CD22", "PAX5", "CD19"
    ),
    "Plasma_B" = c(
        "IGHG1", "MZB1", "SDC1", "JCHAIN",
        "XBP1", "SLAMF7", "IGHM"
    ),

    # --- Other ---
    "Mast_Cell" = c(
        "KIT", "GATA2", "HPGDS", "CPA3",
        "TPSD1", "SLC18A2", "IL1RL1", "CTSG"
    ),
    "Neural_Schwann" = c(
        "MPZ", "SOX10", "GAP43", "PMP22",
        "PLP1", "S100B", "NCAM1", "SOX2"
    )
)

# Compact panel: 1–2 representative genes per lineage
MARKER_COMPACT <- c(
    "TMPRSS2", "NKX3-1",             # Epithelial luminal
    "KRT5", "SOX9",                   # Epithelial basal / club
    "COL1A1", "ACTA2",                # Stromal
    "PECAM1", "VWF",                  # Endothelial
    "CD68", "CD3E",                   # Myeloid / T cell
    "MS4A1", "CD79A",                 # B cell
    "KIT", "MPZ", "SOX10"             # Mast / Neural
)

# ======================== Per-Sample Annotation ========================
rds_files <- list.files(DATA_DIR,
                        pattern = "^purified_xe_output-.*\\.rds$",
                        full.names = TRUE)

if (length(rds_files) == 0) {
    stop("No purified_xe_output-*.rds files found in ", DATA_DIR)
}

message(sprintf("Found %d purified RDS file(s)", length(rds_files)))

for (rds_path in rds_files) {

    # Extract sample ID from filename
    sample_id <- gsub("^purified_xe_output-|\\.rds$", "",
                      basename(rds_path))
    message(sprintf("\n===== %s =====", sample_id))

    obj <- readRDS(rds_path)
    cell_types <- unique(obj$first_type)
    message(sprintf("  Cells: %d   Types: %d   Genes: %d",
                    ncol(obj), length(cell_types), nrow(obj)))

    # ---- 1. UMAP colored by first_type ----
    p1 <- DimPlot(obj, group.by = "first_type", label = TRUE,
                  repel = TRUE, raster = TRUE) +
          theme(aspect.ratio = 1) +
          ggtitle(paste(sample_id, "— SPLIT first_type"))
    ggsave(file.path(OUTPUT_DIR, paste0(sample_id, "_01_umap_first_type.png")),
           plot = p1, width = 10, height = 8, dpi = 150)

    # ---- 2. Compact marker FeaturePlot ----
    genes_present <- intersect(MARKER_COMPACT, rownames(obj))
    if (length(genes_present) > 0) {
        p2 <- FeaturePlot(obj, features = genes_present,
                          raster = FALSE, order = TRUE, ncol = 5) +
              plot_annotation(title = paste(sample_id, "— Marker Genes"))
        n_genes <- length(genes_present)
        ggsave(file.path(OUTPUT_DIR, paste0(sample_id, "_02_markers.png")),
               plot = p2,
               width = min(20, 4 * min(5, n_genes)),
               height = max(8, 4 * ceiling(n_genes / 5)),
               dpi = 150)
    }

    # ---- 3. DotPlot: markers × first_type ----
    all_markers <- unique(unlist(MARKER_PANELS))
    markers_avail <- intersect(all_markers, rownames(obj))
    if (length(markers_avail) > 3) {
        p3 <- DotPlot(obj, features = markers_avail,
                      group.by = "first_type") +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
              ggtitle(paste(sample_id, "— Marker DotPlot"))
        ggsave(file.path(OUTPUT_DIR,
                         paste0(sample_id, "_03_dotplot.png")),
               plot = p3,
               width = max(12, length(markers_avail) * 0.3),
               height = 8, dpi = 150)
    }

    # ---- 4. Cell type proportion table ----
    cell_stats <- as.data.frame(table(obj$first_type))
    colnames(cell_stats) <- c("CellType", "Count")
    cell_stats$Pct <- round(
        cell_stats$Count / sum(cell_stats$Count) * 100, 2)
    cell_stats <- cell_stats[order(-cell_stats$Count), ]
    write.csv(cell_stats,
              file.path(OUTPUT_DIR,
                        paste0(sample_id, "_04_celltype_proportions.csv")),
              row.names = FALSE)

    # Print summary
    message("  Top 5 cell types:")
    for (j in seq_len(min(5, nrow(cell_stats)))) {
        message(sprintf("    %-20s %6d (%5.1f%%)",
                cell_stats$CellType[j],
                cell_stats$Count[j],
                cell_stats$Pct[j]))
    }

    # ---- 5. Epithelial sub-clustering (if Epi detected) ----
    epi_labels <- c("Luminal", "Club", "Basal")
    epi_present <- intersect(epi_labels, cell_types)
    if (length(epi_present) > 0) {
        message(sprintf("  Sub-clustering Epithelial (%s)...",
                        paste(epi_present, collapse=", ")))
        epi <- subset(obj, subset = first_type %in% epi_present)
        epi <- SCTransform(epi, assay = "Xenium") %>%
               RunPCA() %>% RunUMAP(dims = 1:30)

        p5 <- DimPlot(epi, group.by = "first_type", label = TRUE,
                      repel = TRUE) +
              theme(aspect.ratio = 1) +
              ggtitle(paste(sample_id, "— Epithelial Subtypes"))
        ggsave(file.path(OUTPUT_DIR,
                         paste0(sample_id, "_05_epi_subtypes.png")),
               plot = p5, width = 8, height = 7, dpi = 150)
    }
}

message("\n=== Annotation complete ===")
