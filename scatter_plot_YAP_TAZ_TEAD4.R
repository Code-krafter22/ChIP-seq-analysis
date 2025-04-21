#Load Required Libraries
library(GenomicRanges)
library(rtracklayer) # for reading in bed file
library(here)
library(dplyr)
library(ggplot2)

TAZ_peaks <- import(here("C:\\Users\\DELL\\Downloads\\acg_midterm\\new_data\\final_peak\\TAZ_macs3_peak\\TAZmacs3_peaks.narrowPeak"))
YAP_peaks <- import(here("C:\\Users\\DELL\\Downloads\\acg_midterm\\new_data\\final_peak\\YAP_macs3_peak\\YAP_macs3_peaks.narrowPeak"))

common_levels <- intersect(seqlevels(TAZ_peaks), seqlevels(YAP_peaks))

TAZ_peaks <- keepSeqlevels(TAZ_peaks, common_levels, pruning.mode = "coarse")
YAP_peaks <- keepSeqlevels(YAP_peaks, common_levels, pruning.mode = "coarse")

YAP_overlap_TAZ_peaks<- subsetByOverlaps(YAP_peaks, TAZ_peaks)
length(YAP_overlap_TAZ_peaks)

# use rtracklayer to write the GenomicRanges object to file
export(YAP_overlap_TAZ_peaks_overlap_TEAD4, 
       con = here("C:\\Users\\DELL\\Downloads\\acg_midterm\\new_data\\final_peak\\YAP_TAZ_TEAD4_common.bed"))

# after getting the .bed file count the number of reads from bam files with bedtools
# bedtools multicov -bams SRR1810900.bam SRR1810912.bam SRR1810918.bam -bed YAP_TAZ_TEAD4_common.bed > YAP_TAZ_TEAD4_counts.tsv

library(readr)
counts<- read_tsv(here("C:\\Users\\DELL\\Downloads\\acg_midterm\\new_data\\final_peak\\YAP_TAZ_TEAD4_counts.tsv"), col_names = FALSE)
colnames(counts)<- c("chr", "start", "end", "name", "score", "value", "YAP1", "TAZ", "TEAD4")

head(counts)

counts<- counts %>%
  mutate(YAP1 = YAP1/23653961 * 10^6,
         TAZ = TAZ/26789648 * 10^6,
         TEAD4 = TEAD4/34332907 * 10^6)

head(counts)

library(ggpmisc)
ggplot(counts, aes(x = TEAD4, y = YAP1)) +
  geom_point(color = "seagreen") +
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
  ylab("YAP1 signal")

#calculate correlation coefficient for YAP vs TEAD4

correlation_coefficient <- cor(log2(counts$TEAD4), log2(counts$YAP1))
#[1] 0.8095894

R_square <- correlation_coefficient^2
#[1] 0.655435


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

#calculate correlation coefficient for YAP vs TEAD4

correlation_coefficient <- cor(log2(counts$TEAD4), log2(counts$TAZ))
#[1] 0.8179775

R_square <- correlation_coefficient^2
#[1] 0.6690872
