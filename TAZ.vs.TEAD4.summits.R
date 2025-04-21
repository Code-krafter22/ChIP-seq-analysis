library(GenomicRanges)
library(rtracklayer)
library(here)

#Import summit data
TAZ_summit <- import(here("C:\\Users\\DELL\\Downloads\\acg_midterm\\new_data\\final_peak\\TAZ_macs3_peak\\TAZmacs3_summits.bed"))
TAZ_summit <- TAZ_summit[TAZ_summit$name %in% TAZ_overlap_YAP_peaks$name]
TEAD4_summit <- import(here("C:\\Users\\DELL\\Downloads\\acg_midterm\\new_data\\final_peak\\TEAD4_macs3_peak\\TEAD4_macs3_summits.bed"))
TEAD4_summit

#expand the TAZ summit to a 500bp window
TAZ_500bp_window<- resize(TAZ_summit, width = 500, fix="center")

#Finds overlapping peaks and computes distances between TEAD4 and TAZ summits.
hits<- findOverlaps(TEAD4_summit, TAZ_500bp_window)

hits

summit_distance<- distance(TEAD4_summit[queryHits(hits)], TAZ_summit[subjectHits(hits)])

table(summit_distance)

TEAD4_summit[queryHits(hits)][summit_distance ==0]

TAZ_summit[subjectHits(hits)][summit_distance ==0]

# Compute signed distances
signed_distance <- function(A, B) {
  # Compute unsigned distance
  dist <- distance(A, B)
  
  # Determine signs based on whether A precedes or follows B
  sign <- ifelse(start(A) < start(B), -1, 1)
  
  # Apply sign to distance
  dist * sign
}

library(dplyr)
library(ggplot2)
summit_distance<- signed_distance(TEAD4_summit[queryHits(hits)],TAZ_summit[subjectHits(hits)])

distance_df<- table(summit_distance) %>%
  tibble::as_tibble() 

distance_df

distance_df %>%
  mutate(summit_distance = as.numeric(summit_distance)) %>%
  arrange(summit_distance) %>%
  ggplot(aes(x=summit_distance, y = n)) +
  geom_line()

df_binned <- distance_df %>%
  mutate(summit_distance = as.numeric(summit_distance)) %>%
  arrange(summit_distance) %>%
  mutate(bin = floor(summit_distance / 5) * 5) %>%  # Create bins by grouping every 5 bp
  group_by(bin) %>%
  summarise(n = mean(n, na.rm = TRUE))  # Calculate average 'n' for each bin

# View the binned dataframe
print(df_binned)

df_binned %>%
  ggplot(aes(x=bin, y = n)) +
  geom_line() +
  scale_x_continuous(breaks = c(-250, 0, 250)) +
  xlab("distance to the summit \nof TAZ peaks (bp)") +
  ylab("peak density") +
  theme_classic(base_size = 14)
