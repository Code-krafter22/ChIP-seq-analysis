# ==========================
# ChIP-seq Overlap Analysis
# ==========================

# Load Required Libraries
library(GenomicRanges)
library(rtracklayer)    # For importing and exporting BED files
library(here)           # To construct file paths
library(dplyr)          # Data manipulation
library(ggplot2)        # Visualization

# -------------------------
# Import Peak Files
# -------------------------

# Import narrowPeak files for TAZ and YAP from MACS3
TAZ_peaks <- import(here("C:\\Users\\DELL\\Downloads\\acg_midterm\\new_data\\final_peak\\TAZ_macs3_peak\\TAZmacs3_peaks.narrowPeak"))
YAP_peaks <- import(here("C:\\Users\\DELL\\Downloads\\acg_midterm\\new_data\\final_peak\\YAP_macs3_peak\\YAPmacs3_peaks.narrowPeak"))

# Ensure both peak sets are on the same set of chromosomes (seqlevels)
common_levels <- intersect(seqlevels(TAZ_peaks), seqlevels(YAP_peaks))

TAZ_peaks <- keepSeqlevels(TAZ_peaks, common_levels, pruning.mode = "coarse")
YAP_peaks <- keepSeqlevels(YAP_peaks, common_levels, pruning.mode = "coarse")

# -------------------------
# Find Overlapping Peaks
# -------------------------

# Identify YAP peaks that overlap with TAZ peaks
YAP_overlap_TAZ_peaks <- subsetByOverlaps(YAP_peaks, TAZ_peaks)
length(YAP_overlap_TAZ_peaks)  # View number of overlapping peaks

# Export the overlapping peaks for downstream analysis (e.g., TEAD4 intersection)
export(YAP_overlap_TAZ_peaks_overlap_TEAD4, 
       con = here("C:\\Users\\DELL\\Downloads\\acg_midterm\\new_data\\final_peak\\YAP_TAZ_TEAD4_common.bed"))

# BEDTools command to generate read counts:
# bedtools multicov -bams SRR1810900.bam SRR1810912.bam SRR1810918.bam \
# -bed YAP_TAZ_TEAD4_common.bed > YAP_TAZ_TEAD4_counts.tsv

# -------------------------
# Read and Normalize Counts
# -------------------------

library(readr)
# Read in read counts file (no headers)
counts <- read_tsv(here("C:\\Users\\DELL\\Downloads\\acg_midterm\\new_data\\final_peak\\YAP_TAZ_TEAD4_counts.tsv"), 
                   col_names = FALSE)

# Assign appropriate column names
colnames(counts) <- c("chr", "start", "end", "name", "score", "value", "YAP1", "TAZ", "TEAD4")

# Normalize read counts to CPM (counts per million)
counts <- counts %>%
  mutate(
    YAP1 = YAP1 / 23653961 * 1e6,
    TAZ = TAZ / 26789648 * 1e6,
    TEAD4 = TEAD4 / 34332907 * 1e6
  )

# -------------------------
# Plot YAP1 vs TEAD4
# -------------------------

library(ggpmisc)

ggplot(counts, aes(x = TEAD4, y = YAP1)) +
  geom_point(color = "seagreen") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add linear regression line
  stat_poly_eq(  # Show R-squared on the plot
    aes(label = after_stat(rr.label)),
    formula = y ~ x,
    parse = TRUE,
    color = "black"
  ) +
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2') +
  theme_classic(base_size = 14) +
  xlab("TEAD4 signal") +
  ylab("YAP1 signal")

# Calculate and display correlation and R-squared
correlation_coefficient <- cor(log2(counts$TEAD4), log2(counts$YAP1))
R_square <- correlation_coefficient^2
# Print values
correlation_coefficient  # ~0.81
R_square                 # ~0.66

# -------------------------
# Plot TAZ vs TEAD4
# -------------------------

ggplot(counts, aes(x = TEAD4, y = TAZ)) +
  geom_point(color = "tomato") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_poly_eq(
    aes(label = after_stat(rr.label)),
    formula = y ~ x,
    parse = TRUE,
    color = "black"
  ) +
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2') +
  theme_classic(base_size = 14) +
  xlab("TEAD4 signal") +
  ylab("TAZ signal")

# Calculate and display correlation and R-squared
correlation_coefficient <- cor(log2(counts$TEAD4), log2(counts$TAZ))
R_square <- correlation_coefficient^2
# Print values
correlation_coefficient  # ~0.82
R_square                 # ~0.67
