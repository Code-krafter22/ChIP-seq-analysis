# We are trying to plot figures from a research paper : 

# The experiment runs ChiP-seq experiment. 

# So far, we have taken the sequencing data from ENA, and we completed our pre-processing. 
# Ultimately, we have bigwig files, peak files and bed files, which are essential formats of ChIP-Seq data analysis. 
# We will now use these data for visualization.  

library(GenomicRanges)
library(rtracklayer)
library(here)

TAZ_peaks <- import(here("C:\\Users\\DELL\\Downloads\\acg_midterm\\new_data\\final_peak\\TAZ_macs3_peak\\TAZmacs3_peaks.narrowPeak"))
YAP_peaks <- import(here("C:\\Users\\DELL\\Downloads\\acg_midterm\\new_data\\final_peak\\YAP_macs3_peak\\YAP_macs3_peaks.narrowPeak"))

TAZ_peaks
YAP_peaks

seqlevels(TAZ_peaks)
seqlevels(YAP_peaks)

common_levels <- intersect(seqlevels(TAZ_peaks), seqlevels(YAP_peaks))

TAZ_peaks <- keepSeqlevels(TAZ_peaks, common_levels, pruning.mode = "coarse")
YAP_peaks <- keepSeqlevels(YAP_peaks, common_levels, pruning.mode = "coarse")

TAZ_overlap_YAP_peaks <-subsetByOverlaps(TAZ_peaks, YAP_peaks)
length(TAZ_overlap_YAP_peaks)

YAP_overlap_TAZ_peaks <-subsetByOverlaps(YAP_peaks, TAZ_peaks)
length(YAP_overlap_TAZ_peaks)


# install Vennerable from :
#     devtools::install_github("js229/Vennerable")
# It may require you to download additional packages such as:
#     install("RBGL")
library(Vennerable)

n_YAP <- length(YAP_peaks)
n_TAZ <- length(TAZ_peaks)
n_overlap <- length(YAP_overlap_TAZ_peaks)

# Figure 1a :  venn-diagram showing the overlapping peak number of YAP and TAZ 
venn_data <- Venn(SetNames = c("YAP", "TAZ"),
                  Weight = c(
                    "10" = n_YAP - n_overlap , # unique to A n(A) not "n only (A)"
                    "01" = n_TAZ - n_overlap, # unique to B
                    "11" = n_overlap # Intersection
                  ))

plot(venn_data)

# Figure 1b:  venn-diagram showing the overlapping peak number of YAP/TAZ and TEAD4
TEAD4_peaks<- import(here("C:\\Users\\DELL\\Downloads\\acg_midterm\\new_data\\final_peak\\TEAD4_macs3_peak\\TEAD4_macs3_peaks.narrowPeak"))
TEAD4_peaks

YAP_overlap_TAZ_peaks_overlap_TEAD4 <- subsetByOverlaps(YAP_overlap_TAZ_peaks, TEAD4_peaks)

n_YAP_TAZ <- length(YAP_overlap_TAZ_peaks)
n_TEAD4 <- length(TEAD4_peaks)
n_overlap2 <- length(YAP_overlap_TAZ_peaks_overlap_TEAD4)
n_overlap2

venn_data2 <- Venn(SetNames = c("YAP/TAZ", "TEAD4"),
                   Weight = c(
                     "10" = n_YAP_TAZ-n_overlap2,
                     "01" = n_TEAD4-n_overlap2,
                     "11" = n_overlap2
                   ))

plot(venn_data2)
