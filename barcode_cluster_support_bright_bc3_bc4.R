library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)

# ---------------------------------
# paths
# ---------------------------------
dir_bc4 <- "/Users/anissab/Desktop/illumina_shuffled_multi-bin_barcode_analysis/barcodes_and_mutants/batch_barcode_comparison_results_clusterlevel_mix_vs_bc4"

dir_bc3 <- "/Users/anissab/Desktop/illumina_shuffled_multi-bin_barcode_analysis/barcodes_and_mutants/batch_barcode_comparison_results_clusterlevel_mix_vs_bc3"

output_file <- "/Users/anissab/Desktop/illumina_shuffled_multi-bin_barcode_analysis/barcodes_and_mutants/mix_bc3_bc4_barcode_cluster_support.csv"

# ---------------------------------
# helper: find paired pacbio file
# ---------------------------------
find_pacbio_file <- function(file_base, folder_path, suffix_name) {
  
  candidates <- c(
    # case like:
    # matched_clusters_mix_bc4 -> matched_pacbio_mix_bc4
    sub("matched_clusters", "matched_pacbio", file_base),
    
    # case like:
    # matched_clusters.tsv -> matched_pacbio_bc3.tsv
    sub("_matched_clusters\\.tsv$", paste0("_matched_pacbio_", suffix_name, ".tsv"), file_base),
    
    # fallback generic
    sub("_matched_clusters\\.tsv$", "_matched_pacbio.tsv", file_base)
  )
  
  candidates <- unique(candidates)
  
  for (cand in candidates) {
    cand_path <- file.path(folder_path, cand)
    if (file.exists(cand_path)) {
      return(cand)
    }
  }
  
  return(NULL)
}

# ---------------------------------
# helper: process one folder
# ---------------------------------
process_cluster_folder <- function(folder_path, suffix_name) {
  
  cluster_files <- list.files(
    folder_path,
    pattern = "matched_clusters.*\\.tsv$",
    full.names = TRUE
  )
  
  if (length(cluster_files) == 0) {
    stop(paste("No matched_clusters files found in:", folder_path))
  }
  
  message("Found ", length(cluster_files), " matched_clusters files in ", folder_path)
  
  per_file_results <- lapply(cluster_files, function(cluster_file) {
    
    file_base <- basename(cluster_file)
    pacbio_file_base <- find_pacbio_file(file_base, folder_path, suffix_name)
    
    if (is.null(pacbio_file_base)) {
      warning(
        paste0(
          "Missing matched_pacbio file for: ", file_base,
          "\nTried suffix_name = ", suffix_name
        )
      )
      return(NULL)
    }
    
    pacbio_file <- file.path(folder_path, pacbio_file_base)
    
    message("Pairing:")
    message("  clusters: ", file_base)
    message("  pacbio:   ", pacbio_file_base)
    
    clusters <- read_tsv(cluster_file, show_col_types = FALSE)
    pacbio <- read_tsv(pacbio_file, show_col_types = FALSE)
    
    required_cluster_cols <- c("representative_barcode", "matched_pacbio_barcodes")
    missing_cluster_cols <- setdiff(required_cluster_cols, colnames(clusters))
    if (length(missing_cluster_cols) > 0) {
      warning(
        paste0(
          "Skipping ", file_base,
          " because it is missing columns: ",
          paste(missing_cluster_cols, collapse = ", ")
        )
      )
      return(NULL)
    }
    
    required_pacbio_cols <- c("barcode", "aa_seq")
    missing_pacbio_cols <- setdiff(required_pacbio_cols, colnames(pacbio))
    if (length(missing_pacbio_cols) > 0) {
      warning(
        paste0(
          "Skipping paired pacbio file for ", file_base,
          " because it is missing columns: ",
          paste(missing_pacbio_cols, collapse = ", ")
        )
      )
      return(NULL)
    }
    
    cluster_long <- clusters %>%
      filter(!is.na(matched_pacbio_barcodes), matched_pacbio_barcodes != "") %>%
      separate_rows(matched_pacbio_barcodes, sep = ",") %>%
      mutate(
        representative_barcode = trimws(representative_barcode),
        matched_pacbio_barcodes = trimws(matched_pacbio_barcodes)
      ) %>%
      filter(
        !is.na(representative_barcode), representative_barcode != "",
        !is.na(matched_pacbio_barcodes), matched_pacbio_barcodes != ""
      )
    
    pacbio_clean <- pacbio %>%
      select(barcode, aa_seq) %>%
      mutate(
        barcode = trimws(barcode),
        aa_seq = trimws(aa_seq)
      ) %>%
      filter(!is.na(barcode), barcode != "", !is.na(aa_seq), aa_seq != "") %>%
      distinct(barcode, aa_seq)
    
    joined <- cluster_long %>%
      left_join(pacbio_clean, by = c("matched_pacbio_barcodes" = "barcode")) %>%
      filter(!is.na(aa_seq), aa_seq != "")
    
    if (nrow(joined) == 0) {
      warning("No joined rows for: ", file_base)
      return(NULL)
    }
    
    joined %>%
      group_by(aa_seq) %>%
      summarise(
        n_barcode_clusters = n_distinct(representative_barcode),
        .groups = "drop"
      ) %>%
      mutate(source_file = file_base)
  })
  
  per_file_results <- per_file_results[!sapply(per_file_results, is.null)]
  
  if (length(per_file_results) == 0) {
    stop(paste("No usable results generated for folder:", folder_path))
  }
  
  combined <- bind_rows(per_file_results)
  
  folder_summary <- combined %>%
    group_by(aa_seq) %>%
    summarise(
      n_barcode_clusters = sum(n_barcode_clusters, na.rm = TRUE),
      .groups = "drop"
    )
  
  colnames(folder_summary)[2] <- paste0("n_barcode_clusters_", suffix_name)
  
  return(folder_summary)
}

# ---------------------------------
# run bc4 and bc3
# ---------------------------------
bc4_support <- process_cluster_folder(dir_bc4, "bc4")
bc3_support <- process_cluster_folder(dir_bc3, "bc3")

# ---------------------------------
# merge and total
# ---------------------------------
support_merged <- full_join(bc4_support, bc3_support, by = "aa_seq") %>%
  mutate(
    n_barcode_clusters_bc4 = replace_na(n_barcode_clusters_bc4, 0),
    n_barcode_clusters_bc3 = replace_na(n_barcode_clusters_bc3, 0),
    n_barcode_clusters_total = n_barcode_clusters_bc4 + n_barcode_clusters_bc3
  ) %>%
  arrange(desc(n_barcode_clusters_total))

# ---------------------------------
# save
# ---------------------------------
write_csv(support_merged, output_file)

cat("Done.\n")
cat("Wrote:", output_file, "\n\n")

cat("Summary of total barcode clusters:\n")
print(summary(support_merged$n_barcode_mix_clusters_total))

cat("\nTop supported sequences:\n")
print(head(support_merged, 10))

cat("\nNumber of sequences in merged support table:\n")
print(nrow(support_merged))
