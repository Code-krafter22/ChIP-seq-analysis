library(rtracklayer) # for reading in bed file
library(here)
library(dplyr)
library(ggplot2)

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


H3K4me1<- import(here("data/public_data/H3K4me1.bed"))
H3K4me3<- import(here("data/public_data/H3K4me3.bed"))
H3K27ac<- import(here("data/public_data/H3K27ac.bed"))

active_enhancers<- subsetByOverlaps(H3K4me1, H3K27ac)
inactive_enhancers<- subsetByOverlaps(H3K4me1, H3K27ac, invert=TRUE)

promoters<- subsetByOverlaps(H3K4me3, H3K4me1, invert=TRUE)

n_active_enhancers<- subsetByOverlaps(YAP_overlap_TAZ_peaks_overlap_TEAD4,
                                      active_enhancers) %>%
  length()

n_inactive_enhancers<- subsetByOverlaps(YAP_overlap_TAZ_peaks_overlap_TEAD4,
                                        inactive_enhancers) %>%
  length()

n_promoters<- subsetByOverlaps(YAP_overlap_TAZ_peaks_overlap_TEAD4, 
                               promoters) %>%
  length()

n_unclassified<- length(YAP_overlap_TAZ_peaks_overlap_TEAD4) - n_active_enhancers -
  n_inactive_enhancers - n_promoters

annotation_df<- data.frame(category = c("active_enhancers", "inactive_enhancers",
                                        "promoters", "unclassified"),
                           peak_number = c(n_active_enhancers, n_inactive_enhancers, 
                                           n_promoters, n_unclassified))


annotation_df

library(ggplot2)

ggplot(annotation_df, aes(x = "", y = peak_number, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() + # Remove unnecessary axes
  labs(title = "YAP/TAZ/TEAD4 peaks") +
  scale_fill_brewer(palette = "Set3") 