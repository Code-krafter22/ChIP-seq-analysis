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

YAP_summit<- import(here("C:\\Users\\DELL\\Downloads\\acg_midterm\\new_data\\final_peak\\YAP_macs3_peak\\YAP_macs3_summits.bed"))
YAP_summit

H3K4me1<- import(here("data/public_data/H3K4me1.bed"))
H3K4me3<- import(here("data/public_data/H3K4me3.bed"))
H3K27ac<- import(here("data/public_data/H3K27ac.bed"))

enhancers<- subsetByOverlaps(H3K4me1, H3K4me3, invert=TRUE)

promoters<- subsetByOverlaps(H3K4me3, H3K4me1, invert=TRUE)

  
YAP1_enhancers<- subsetByOverlaps(YAP_overlap_TAZ_peaks_overlap_TEAD4, enhancers) 

YAP1_promoters<- subsetByOverlaps(YAP_overlap_TAZ_peaks_overlap_TEAD4, promoters) 
YAP1_enhancers$name %>% head()

YAP_summit_enhancer<- YAP_summit[YAP_summit$name %in% YAP1_enhancers$name]
YAP_summit_promoter<- YAP_summit[YAP_summit$name %in% YAP1_promoters$name]

# combine them
anchors<- c(YAP_summit_promoter, YAP_summit_enhancer) 


YAP1_bw<- import(here("data/fastq/YAP.bw"))
TAZ_bw<- import(here("data/fastq/TAZ.bw"))
TEAD4_bw<- import(here("data/fastq/TEAD4.bw"))

YAP1_bw

# BiocManager::install("EnrichedHeatmap")
library(EnrichedHeatmap)
# extend 1000 bp on each side and use 50bp bin
mat1<- normalizeToMatrix(YAP1_bw, anchors, value_column = "score",
                         extend= 1000, mean_mode = "w0", w=50)

mat2<- normalizeToMatrix(TAZ_bw, anchors, value_column = "score",
                         extend= 1000, mean_mode = "w0", w=50)

mat3<- normalizeToMatrix(TEAD4_bw, anchors, value_column = "score",
                         extend= 1000, mean_mode = "w0", w=50)

dim(mat1)

dim(mat2)

dim(mat3)

mat1[1:5, 1:40 ]

quantile(mat1, c(0.1,0.25,0.5,0.9,1))

quantile(mat2, c(0.1,0.25,0.5,0.9,1))

quantile(mat3, c(0.1,0.25,0.5,0.9,1))

col_fun<- circlize::colorRamp2(c(0, 20), c("white", "red"))

partition<- c(rep("promoters", length(YAP1_promoters)),
              rep("enhancers", length(YAP1_enhancers)))

# change the factor level so promoters come first
partition<- factor(partition, levels=c("promoters", "enhancers"))

partition_hp<- Heatmap(partition, col=structure(2:3, names = c("promoters", "enhancers")), 
                       name = "partition",
                       show_row_names = FALSE, width=unit(3,'mm'))

partition_hp

draw(ht_list, split= partition, main_heatmap =2)

