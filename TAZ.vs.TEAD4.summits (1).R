
TAZ_summit <- import(here("~/Desktop/Computational_Genomics/Midterm Presentation/TAZ_peak/TAZ_summits.bed"))
TAZ_summit <- TAZ_summit[TAZ_summit$name %in% TAZ_overlap_YAP_peaks$name]
TEAD4_summit <- import(here("~/Desktop/Computational_Genomics/Midterm Presentation/TEAD4_peak/TEAD4_summits.bed"))
TEAD4_summit


#expand the TAZ summit to a 500bp window
TAZ_500bp_window<- resize(TAZ_summit, width = 500, fix="center")

hits<- findOverlaps(TEAD4_summit, TAZ_500bp_window)

# a hits object with the indices of the overlapping query and subject
hits


summit_distance<- distance(TEAD4_summit[queryHits(hits)], TAZ_summit[subjectHits(hits)])

table(summit_distance)

TEAD4_summit[queryHits(hits)][summit_distance ==0]

TAZ_summit[subjectHits(hits)][summit_distance ==0]



#use revised distance() function to return negative values 
#when TEAD4 summit precede the TAZ summit and positive values when TEAD4 summit follows TAZ summit.

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
summit_distance<- signed_distance(TEAD4_summit[queryHits(hits)],
                                  TAZ_summit[subjectHits(hits)])

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


