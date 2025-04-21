# ===============================================
# TSS Proximity Analysis for YAP, TAZ, and TEAD4
# ===============================================
# This script calculates the distance of ChIP-seq peaks (YAP, TAZ, TEAD4)
# to the nearest transcription start sites (TSS), categorizes them by
# distance range, and visualizes their distributions.

# ------------------------------------------------
# Load required libraries and genome annotations
# ------------------------------------------------

# Install if needed:
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # Transcript annotation for hg38
library(GenomicRanges)                      # Working with genomic coordinates
library(GenomicFeatures)                    # TSS/promoter extraction
library(dplyr)                              # Data manipulation
library(ggplot2)                            # Plotting

# ------------------------------------------------
# Get the transcript start sites (TSS)
# ------------------------------------------------

# Retrieve all transcript annotations from hg38
hg38_transcripts <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Define TSS as 1bp region at the start of each transcript
tss_gr <- promoters(hg38_transcripts, upstream = 0, downstream = 1)

# ------------------------------------------------
# Calculate distance from peaks to nearest TSS
# ------------------------------------------------

# Calculate distance for each peak to its nearest TSS
distance_to_tss <- distanceToNearest(YAP_peaks, tss_gr)

# View result object and extract metadata
distance_to_tss
mcols(distance_to_tss)
head(mcols(distance_to_tss)$distance)

# Extract distance values for each factor
YAP_dist <- mcols(distanceToNearest(YAP_peaks, tss_gr))$distance
TAZ_dist <- mcols(distanceToNearest(TAZ_peaks, tss_gr))$distance
TEAD4_dist <- mcols(distanceToNearest(TEAD4_peak, tss_gr))$distance

# ------------------------------------------------
# Combine distances into a single data frame
# ------------------------------------------------

tss_distance_df <- bind_rows(
  data.frame(factor = "YAP", distance = YAP_dist),
  data.frame(factor = "TAZ", distance = TAZ_dist),
  data.frame(factor = "TEAD4", distance = TEAD4_dist)
)

# Preview: categorize distances into bins
tss_distance_df %>%
  mutate(category = case_when(
    distance < 1000 ~ "<1kb",
    distance >= 1000 & distance < 10000 ~ "1-10kb",
    distance >= 10000 & distance <= 100000 ~ "10-100kb",
    distance > 100000 ~ ">100kb"
  )) %>%
  head()

# ------------------------------------------------
# Count peaks per category and factor
# ------------------------------------------------

counts_per_category <- tss_distance_df %>%
  mutate(category = case_when(
    distance < 1000 ~ "<1kb",
    distance >= 1000 & distance < 10000 ~ "1-10kb",
    distance >= 10000 & distance <= 100000 ~ "10-100kb",
    distance > 100000 ~ ">100kb"
  )) %>%
  group_by(factor, category) %>%
  count()

counts_per_category  # View counts by group and distance category

# ------------------------------------------------
# Total number of peaks per factor (for percentage calc)
# ------------------------------------------------

total_counts <- tss_distance_df %>%
  mutate(category = case_when(
    distance < 1000 ~ "<1kb",
    distance >= 1000 & distance < 10000 ~ "1-10kb",
    distance >= 10000 & distance <= 100000 ~ "10-100kb",
    distance > 100000 ~ ">100kb"
  )) %>%
  count(factor, name = "total")

total_counts

# ------------------------------------------------
# Merge and visualize
# ------------------------------------------------

# Ensure factor levels are ordered for consistent plotting
merged_df$category <- factor(merged_df$category, 
                             levels = c("<1kb", "1-10kb", "10-100kb", ">100kb"))

# Create stacked bar plot of distance categories by factor
merged_df %>%
  mutate(Percentage = n / total * 100) %>%
  ggplot(aes(x = factor, y = Percentage, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Distance to TSS",
    x = "Group",
    y = "Percentage"
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF")) +  # Custom blue shades
  