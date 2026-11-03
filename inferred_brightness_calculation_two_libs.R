library(dplyr)
library(readr)
library(tidyr)

setwd("/Users/anissab/Desktop/illumina_shuffled_multi-bin_barcode_analysis/fluorescence_normalization")

# ----------------------------
# 1. File paths
# ----------------------------
base_dir <- "/Users/anissab/Desktop/illumina_shuffled_multi-bin_barcode_analysis/fluorescence_normalization"
counts_dir <- file.path(base_dir, "mappings_and_counts_bright")

# Sorted brighter-library counts (used for numerator)
sorted_bc3_file <- file.path(
  counts_dir,
  "pacbio_matched_variant_clustered_counts_by_bin_brighter_library_vs_bc3",
  "batch_barcode_comparison_results_clusterlevel_brighter_library_vs_bc3.csv"
)

sorted_bc4_file <- file.path(
  counts_dir,
  "pacbio_matched_variant_clustered_counts_by_bin_brighter_library_vs_bc4",
  "batch_barcode_comparison_results_clusterlevel_brighter_library_vs_bc4.csv"
)

# Presort brighter-library counts
presort_bright_bc3_file <- file.path(
  counts_dir,
  "pacbio_matched_variant_clustered_counts_by_bin_original_brighter_library_vs_bc3",
  "batch_barcode_comparison_results_clusterlevel_original_brighter_library_vs_bc3.csv"
)

presort_bright_bc4_file <- file.path(
  counts_dir,
  "pacbio_matched_variant_clustered_counts_by_bin_original_brighter_library_vs_bc4",
  "batch_barcode_comparison_results_original_brighter_library_vs_bc4.csv"
)

# Presort medium-library counts
presort_medium_bc3_file <- file.path(
  counts_dir,
  "pacbio_matched_variant_clustered_counts_by_bin_original_medium_library_vs_bc3",
  "batch_barcode_comparison_results_original_medium_library_vs_bc3.csv"
)

presort_medium_bc4_file <- file.path(
  counts_dir,
  "pacbio_matched_variant_clustered_counts_by_bin_original_medium_library_vs_bc4",
  "batch_barcode_comparison_results_original_medium_library_vs_bc4.csv"
)

# pooled Medium+Brighter calibration table from FlowCal/raw-bin workflow
bin_mefl_file <- file.path(
  base_dir,
  "Samples_Brighter_plus_Medium_combined_bin_stats_mefl_from_rawbins.csv"
)

# ----------------------------
# 2. Column names
# ----------------------------
bin_cols <- c(
  "B12_S11", "B13_S12", "B14_S13", "B15_S14", "B16_S15",
  "B17_S16", "B18_S17", "B19_S18", "B20_S19"
)

presort_bright_col <- "Brighter-Library_S20"
presort_medium_col <- "Medium-Library_S10"

# ----------------------------
# 3. Sanity checks
# ----------------------------
print(sorted_bc3_file)
print(sorted_bc4_file)
print(presort_bright_bc3_file)
print(presort_bright_bc4_file)
print(presort_medium_bc3_file)
print(presort_medium_bc4_file)
print(bin_mefl_file)

print(file.exists(sorted_bc3_file))
print(file.exists(sorted_bc4_file))
print(file.exists(presort_bright_bc3_file))
print(file.exists(presort_bright_bc4_file))
print(file.exists(presort_medium_bc3_file))
print(file.exists(presort_medium_bc4_file))
print(file.exists(bin_mefl_file))

# ----------------------------
# 4. Read sorted count tables and combine bc3 + bc4
# ----------------------------
sorted_bc3 <- read_csv(sorted_bc3_file, show_col_types = FALSE)
sorted_bc4 <- read_csv(sorted_bc4_file, show_col_types = FALSE)

sorted_bc3 <- sorted_bc3 %>% filter(!is.na(aa_seq), trimws(aa_seq) != "")
sorted_bc4 <- sorted_bc4 %>% filter(!is.na(aa_seq), trimws(aa_seq) != "")

sorted_counts <- bind_rows(sorted_bc3, sorted_bc4) %>%
  group_by(aa_seq) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

missing_bin_cols <- setdiff(bin_cols, colnames(sorted_counts))
if (length(missing_bin_cols) > 0) {
  stop(
    paste(
      "These sorted bin columns are missing from sorted_counts:",
      paste(missing_bin_cols, collapse = ", ")
    )
  )
}

# ----------------------------
# 5. Read presort brighter-library counts
# ----------------------------
presort_bright_bc3 <- read_csv(presort_bright_bc3_file, show_col_types = FALSE)
presort_bright_bc4 <- read_csv(presort_bright_bc4_file, show_col_types = FALSE)

presort_bright_bc3 <- presort_bright_bc3 %>% filter(!is.na(aa_seq), trimws(aa_seq) != "")
presort_bright_bc4 <- presort_bright_bc4 %>% filter(!is.na(aa_seq), trimws(aa_seq) != "")

presort_bright <- bind_rows(presort_bright_bc3, presort_bright_bc4) %>%
  group_by(aa_seq) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

if (!(presort_bright_col %in% colnames(presort_bright))) {
  stop(
    paste(
      "Brighter presort column not found:", presort_bright_col,
      "\nAvailable columns are:\n",
      paste(colnames(presort_bright), collapse = ", ")
    )
  )
}

# ----------------------------
# 6. Read presort medium-library counts
# ----------------------------
presort_medium_bc3 <- read_csv(presort_medium_bc3_file, show_col_types = FALSE)
presort_medium_bc4 <- read_csv(presort_medium_bc4_file, show_col_types = FALSE)

presort_medium_bc3 <- presort_medium_bc3 %>% filter(!is.na(aa_seq), trimws(aa_seq) != "")
presort_medium_bc4 <- presort_medium_bc4 %>% filter(!is.na(aa_seq), trimws(aa_seq) != "")

presort_medium <- bind_rows(presort_medium_bc3, presort_medium_bc4) %>%
  group_by(aa_seq) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

if (!(presort_medium_col %in% colnames(presort_medium))) {
  stop(
    paste(
      "Medium presort column not found:", presort_medium_col,
      "\nAvailable columns are:\n",
      paste(colnames(presort_medium), collapse = ", ")
    )
  )
}

# ----------------------------
# 7. Build 50:50 presort baseline
#    Assumption: assembled library = half brighter + half medium
# ----------------------------
presort_bright_only <- presort_bright %>%
  select(aa_seq, all_of(presort_bright_col)) %>%
  rename(brighter_presort_reads = all_of(presort_bright_col))

presort_medium_only <- presort_medium %>%
  select(aa_seq, all_of(presort_medium_col)) %>%
  rename(medium_presort_reads = all_of(presort_medium_col))

presort_counts <- full_join(presort_bright_only, presort_medium_only, by = "aa_seq") %>%
  mutate(
    brighter_presort_reads = replace_na(brighter_presort_reads, 0),
    medium_presort_reads   = replace_na(medium_presort_reads, 0)
  )

bright_total <- sum(presort_counts$brighter_presort_reads, na.rm = TRUE)
medium_total <- sum(presort_counts$medium_presort_reads, na.rm = TRUE)

presort_counts <- presort_counts %>%
  mutate(
    brighter_presort_prop = ifelse(bright_total > 0, brighter_presort_reads / bright_total, 0),
    medium_presort_prop   = ifelse(medium_total > 0, medium_presort_reads / medium_total, 0),
    presort_norm          = 0.5 * brighter_presort_prop + 0.5 * medium_presort_prop,
    presort_reads         = brighter_presort_reads + medium_presort_reads
  )

# ----------------------------
# 8. Merge sorted + presort baseline
# ----------------------------
flowseq_df <- full_join(sorted_counts, presort_counts, by = "aa_seq") %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

# ----------------------------
# 9. Sorted event totals for bins 12-20
#    Replace with actual brighter-library sorted event totals if needed
# ----------------------------
sorted_events <- tibble(
  bin_col = bin_cols,
  bin = 12:20,
  sorted_events = c(4510577, 22576, 2489, 998, 1638, 2979, 3219, 2289, 1160)
)

# ----------------------------
# 10. Read calibrated bin medians from pooled file
# ----------------------------
bin_mefl <- read_csv(bin_mefl_file, show_col_types = FALSE) %>%
  mutate(
    bin = as.integer(gsub("^P", "", bin)),
    cal_median = median_brightness
  ) %>%
  select(bin, cal_median) %>%
  filter(bin %in% 12:20) %>%
  arrange(bin)

print(bin_mefl)

if (nrow(bin_mefl) != 9) {
  stop("Expected exactly 9 calibrated bins (12-20) in pooled calibration file.")
}

# ----------------------------
# 11. Normalize sorted counts by sequencing depth
# ----------------------------
bin_totals_seq <- flowseq_df %>%
  summarise(across(all_of(bin_cols), ~ sum(.x, na.rm = TRUE)))

norm_df <- flowseq_df %>%
  mutate(across(
    all_of(bin_cols),
    ~ .x / bin_totals_seq[[cur_column()]],
    .names = "{.col}_norm"
  ))

# ----------------------------
# 12. Adjust each bin by fraction of sorted events collected
# ----------------------------
total_sorted_events <- sum(sorted_events$sorted_events, na.rm = TRUE)

event_frac <- sorted_events %>%
  mutate(event_frac = sorted_events / total_sorted_events)

for (i in seq_along(bin_cols)) {
  col_in <- paste0(bin_cols[i], "_norm")
  col_out <- paste0(bin_cols[i], "_adj")
  frac_i <- event_frac$event_frac[event_frac$bin_col == bin_cols[i]]
  
  norm_df[[col_out]] <- norm_df[[col_in]] * frac_i
}

adj_cols <- paste0(bin_cols, "_adj")

# ----------------------------
# 13. Avoid divide-by-zero in presort
# ----------------------------
norm_df <- norm_df %>%
  mutate(presort_norm = ifelse(is.na(presort_norm) | presort_norm == 0, 1e-10, presort_norm))

# ----------------------------
# 14. Fold-change relative to presort
# ----------------------------
for (col in adj_cols) {
  fc_col <- sub("_adj$", "_fc", col)
  norm_df[[fc_col]] <- norm_df[[col]] / norm_df$presort_norm
}

fc_cols <- sub("_adj$", "_fc", adj_cols)

# ----------------------------
# 15. Convert fold-change profile into probability across bins
# ----------------------------
norm_df <- norm_df %>%
  mutate(rowsum_fc = rowSums(across(all_of(fc_cols)), na.rm = TRUE)) %>%
  filter(rowsum_fc > 0)

for (col in fc_cols) {
  a_col <- sub("_fc$", "_A", col)
  norm_df[[a_col]] <- norm_df[[col]] / norm_df$rowsum_fc
}

a_cols <- sub("_fc$", "_A", fc_cols)

# ----------------------------
# 16. Infer fractional median bin
# ----------------------------
bin_labels <- sorted_events$bin

infer_median_bin <- function(probs, bin_labels) {
  probs <- as.numeric(probs)
  
  if (any(is.na(probs))) return(NA_real_)
  
  cs <- cumsum(probs)
  
  if (all(cs < 0.5)) return(max(bin_labels))
  if (all(cs >= 0.5)) return(min(bin_labels))
  
  upper_idx <- which(cs >= 0.5)[1]
  lower_idx <- upper_idx - 1
  
  if (lower_idx < 1) return(bin_labels[upper_idx])
  
  lower_bin <- bin_labels[lower_idx]
  upper_bin <- bin_labels[upper_idx]
  lower_cs <- cs[lower_idx]
  upper_cs <- cs[upper_idx]
  
  frac <- (0.5 - lower_cs) / (upper_cs - lower_cs)
  lower_bin + frac * (upper_bin - lower_bin)
}

norm_df$median_bin <- apply(
  norm_df[, a_cols],
  1,
  function(x) infer_median_bin(as.numeric(x), bin_labels)
)

# ----------------------------
# 17. Convert inferred bin -> MEFL
# ----------------------------
bin_to_mefl <- approxfun(
  x = bin_mefl$bin,
  y = bin_mefl$cal_median,
  rule = 2
)

norm_df$median_MEFL <- bin_to_mefl(norm_df$median_bin)

# ----------------------------
# 18. Add summary columns
# ----------------------------
norm_df <- norm_df %>%
  mutate(
    total_sorted_reads = rowSums(across(all_of(bin_cols)), na.rm = TRUE)
  ) %>%
  select(
    aa_seq,
    all_of(bin_cols),
    brighter_presort_reads,
    medium_presort_reads,
    presort_reads,
    brighter_presort_prop,
    medium_presort_prop,
    presort_norm,
    total_sorted_reads,
    median_bin,
    median_MEFL,
    everything()
  )

# ----------------------------
# 19. Save
# ----------------------------
write_csv(
  norm_df,
  file.path(base_dir, "brighter_flowseq_inferred_brightness_updated.csv")
)

cat("Done. Wrote: brighter_flowseq_inferred_brightness_updated.csv\n")

