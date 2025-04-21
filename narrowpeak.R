# =====================================================
# ChIP-seq Peak Overlap and Venn Diagram Visualization
# =====================================================

# Load Required Libraries
library(GenomicRanges)   # For handling genomic ranges (like peaks)
library(rtracklayer)     # For importing BED/narrowPeak files
library(here)            # For robust file path handling
library(usethis)         # Optional: for project setup and GitHub integration
library(devtools)        # For installing GitHub packages

# -------------------------------
# Load MACS3 narrowPeak files
# -------------------------------

# Import peak files for TAZ and YAP
TAZ_peaks <- import(here("C:\\Users\\DELL\\Downloads\\acg_midterm\\new_data\\final_peak\\TAZ_macs3_peak\\TAZmacs3_peaks.narrowPeak"))
YAP_peaks <- import(here("C:\\Users\\DELL\\Downloads\\acg_midterm\\new_data\\final_peak\\YAP_macs3_peak\\YAPmacs3_peaks.narrowPeak"))

# View peak objects
TAZ_peaks
YAP_peaks

# Optional: view chromosome identifiers in each peak set
# seqlevels(TAZ_peaks)
# seqlevels(YAP_peaks)

# -------------------------------
# Harmonize seqlevels (chromosomes)
# -------------------------------

# Find common chromosomes between TAZ and YAP peaks
common_levels <- intersect(seqlevels(TAZ_peaks), seqlevels(YAP_peaks))

# Keep only the common chromosomes in each peak set
TAZ_peaks <- keepSeqlevels(TAZ_peaks, common_levels, pruning.mode = "coarse")
YAP_peaks <- keepSeqlevels(YAP_peaks, common_levels, pruning.mode = "coarse")

# -------------------------------
# Find overlaps between peak sets
# -------------------------------

# TAZ peaks overlapping with YAP
TAZ_overlap_YAP_peaks <- subsetByOverlaps(TAZ_peaks, YAP_peaks)
length(TAZ_overlap_YAP_peaks)

# YAP peaks overlapping with TAZ (same as above, but from the YAP perspective)
YAP_overlap_TAZ_peaks <- subsetByOverlaps(YAP_peaks, TAZ_peaks)
length(YAP_overlap_TAZ_peaks)

# -------------------------------
# Install & Load Vennerable for Venn diagrams
# -------------------------------

# If not installed, you can install it using:
# install_github("js229/Vennerable")
library(Vennerable)

# Define counts for Venn diagram
n_YAP <- length(YAP_peaks)         # Total YAP peaks
n_TAZ <- length(TAZ_peaks)         # Total TAZ peaks
n_overlap <- length(YAP_overlap_TAZ_peaks)  # Overlapping peaks

# Create Venn diagram object for YAP vs TAZ
venn_data <- Venn(SetNames = c("YAP", "TAZ"),
                  Weight = c(
                    "10" = n_YAP,     # Only in YAP
                    "01" = n_TAZ,     # Only in TAZ
                    "11" = n_overlap  # In both
                  ))

# Plot the Venn diagram
plot(venn_data)

# -------------------------------
# Add TEAD4 peak overlap analysis
# -------------------------------

# Import TEAD4 peak file
TEAD4_peak <- import(here("C:\\Users\\DELL\\Downloads\\acg_midterm\\new_data\\final_peak\\TEAD4_macs3_peak\\TEAD4_macs3_peaks.narrowPeak"))

# Find peaks common across YAP, TAZ, and TEAD4
YAP_overlap_TAZ_peaks_overlap_TEAD4 <- subsetByOverlaps(YAP_overlap_TAZ_peaks, TEAD4_peak)

# Define counts for Venn diagram
n_YAP_TAZ <- length(YAP_overlap_TAZ_peaks)           # Peaks shared by YAP and TAZ
n_TEAD4 <- length(TEAD4_peak)                         # Total TEAD4 peaks
n_overlap2 <- length(YAP_overlap_TAZ_peaks_overlap_TEAD4)  # Overlap with TEAD4

# Create Venn diagram object for (YAP+TAZ) vs TEAD4
venn_data2 <- Venn(SetNames = c("YAP/TAZ", "TEAD4"),
                   Weight = c(
                     "10" = n_YAP_TAZ,    # Only in YAP/TAZ
                     "01" = n_TEAD4,      # Only in TEAD4
                     "11" = n_overlap2    # Shared
                   ))

# Plot the Venn diagram
plot(venn_data2)
