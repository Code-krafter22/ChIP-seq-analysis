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