setwd("/Users/anissab/Desktop/illumina_shuffled_multi-bin_barcode_analysis/barcodes_and_mutants")

# folders that contain the *_matched_pacbio.tsv files you want to combine
folders <- c(
  "batch_barcode_comparison_results_clusterlevel_mix_vs_bc3",
  "batch_barcode_comparison_results_clusterlevel_mix_vs_bc2",
  "batch_barcode_comparison_results_clusterlevel_bright_vs_bc3",
  "batch_barcode_comparison_results_clusterlevel_bright_vs_bc2"
)

# collect all matched_pacbio files from those folders
matched_files <- unlist(lapply(folders, function(folder) {
  list.files(folder, pattern = "_matched_pacbio\\.tsv$", full.names = TRUE)
}))

if (length(matched_files) == 0) {
  stop("No *_matched_pacbio.tsv files found in the listed folders.")
}

# read all aa_seq values
all_seqs <- unlist(lapply(matched_files, function(f) {
  df <- read.delim(f, sep = "\t", stringsAsFactors = FALSE, quote = "")
  if (!"aa_seq" %in% colnames(df)) {
    stop(paste("No aa_seq column found in:", f))
  }
  trimws(df$aa_seq)
}), use.names = FALSE)

# keep non-empty unique sequences only
all_seqs <- all_seqs[!is.na(all_seqs) & all_seqs != ""]
n_unique_genes <- length(unique(all_seqs))

cat("Total unique mapped PacBio sequences across all files:", n_unique_genes, "\n")