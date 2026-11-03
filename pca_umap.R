library(dplyr)
library(ggplot2)

pca_df2 <- pca_df %>%
  mutate(
    bright_top10 = original_log10_MEFL >= quantile(original_log10_MEFL, 0.90, na.rm = TRUE),
    smoothed_top10 = smoothed_log10_MEFL >= quantile(smoothed_log10_MEFL, 0.90, na.rm = TRUE),
    high_support = bright_total > 5
  )

p_top10_orig <- ggplot(pca_df2, aes(PC1, PC2, color = bright_top10)) +
  geom_point(alpha = 0.6, size = 1) +
  theme_classic() +
  labs(
    title = "PCA of ESM2 embeddings",
    subtitle = "Top 10% original brightness highlighted",
    color = "Top 10%"
  )

p_top10_smooth <- ggplot(pca_df2, aes(PC1, PC2, color = smoothed_top10)) +
  geom_point(alpha = 0.6, size = 1) +
  theme_classic() +
  labs(
    title = "PCA of ESM2 embeddings",
    subtitle = "Top 10% smoothed brightness highlighted",
    color = "Top 10%"
  )

p_support_bin <- ggplot(pca_df2, aes(PC1, PC2, color = high_support)) +
  geom_point(alpha = 0.6, size = 1) +
  theme_classic() +
  labs(
    title = "PCA of ESM2 embeddings",
    subtitle = "Support > 5 highlighted",
    color = "Support > 5"
  )

print(p_top10_orig)
print(p_top10_smooth)
print(p_support_bin)

ggplot(pca_df, aes(PC1, PC2, color = aa_length.x)) +
  geom_point(alpha = 0.6, size = 1) +
  theme_classic() +
  labs(
    title = "PCA of ESM2 embeddings",
    subtitle = "Colored by aa_length",
    color = "aa_length.x"
  )
####

library(readr)
library(dplyr)
library(ggplot2)

# ==============================
# File paths
# ==============================
setwd("/Users/anissab/Desktop/illumina_shuffled_multi-bin_barcode_analysis/ESM2_embeddings")

metadata_file   <- "esm2_embeddings_len150_metadata.csv"
brightness_file <- "compare_df_global_bright_prior_smoothed.csv"
embedding_file  <- "esm2_embeddings_len150_matrix.csv"

# ==============================
# Install/load UMAP package
# ==============================
if (!requireNamespace("uwot", quietly = TRUE)) {
  install.packages("uwot")
}
library(uwot)

# ==============================
# Load data
# ==============================
meta <- read_csv(metadata_file, show_col_types = FALSE)
compare_df <- read_csv(brightness_file, show_col_types = FALSE)
embeddings <- as.matrix(read_csv(embedding_file, show_col_types = FALSE))

cat("Embeddings shape:", dim(embeddings), "\n")
cat("Metadata rows:", nrow(meta), "\n")
cat("Brightness rows:", nrow(compare_df), "\n")

# ==============================
# Merge metadata with labels
# ==============================
plot_df <- meta %>%
  rename(aa_length_meta = aa_length) %>%
  inner_join(
    compare_df %>%
      select(aa_seq, aa_length, original_log10_MEFL, smoothed_log10_MEFL, bright_total) %>%
      rename(aa_length_brightness = aa_length),
    by = "aa_seq"
  )

cat("Merged rows:", nrow(plot_df), "\n")

# ==============================
# Align embeddings to merged rows
# ==============================
plot_df <- plot_df %>%
  filter(!is.na(embedding_row)) %>%
  arrange(embedding_row)

embedding_idx <- plot_df$embedding_row + 1   # embedding_row is 0-based
X <- embeddings[embedding_idx, ]

cat("Aligned embedding matrix:", dim(X), "\n")
cat("Aligned plotting rows:", nrow(plot_df), "\n")

# ==============================
# Run UMAP
# ==============================
set.seed(123)

umap_coords <- umap(
  X,
  n_neighbors = 30,
  min_dist = 0.2,
  metric = "euclidean",
  verbose = TRUE
)

umap_df <- plot_df %>%
  mutate(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2]
  )

# ==============================
# Continuous-color UMAP plots
# ==============================
u1 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = original_log10_MEFL)) +
  geom_point(alpha = 0.6, size = 1) +
  theme_classic() +
  labs(
    title = "UMAP of ESM2 embeddings",
    subtitle = "Colored by original_log10_MEFL",
    color = "Original\nlog10 MEFL"
  )

u2 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = smoothed_log10_MEFL)) +
  geom_point(alpha = 0.6, size = 1) +
  theme_classic() +
  labs(
    title = "UMAP of ESM2 embeddings",
    subtitle = "Colored by smoothed_log10_MEFL",
    color = "Smoothed\nlog10 MEFL"
  )

u3 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = log10(bright_total + 1))) +
  geom_point(alpha = 0.6, size = 1) +
  theme_classic() +
  labs(
    title = "UMAP of ESM2 embeddings",
    subtitle = "Colored by support",
    color = "log10(support+1)"
  )

print(u1)
print(u2)
print(u3)

# ==============================
# Binary-highlight UMAP plots
# ==============================
umap_df2 <- umap_df %>%
  mutate(
    bright_top10 = original_log10_MEFL >= quantile(original_log10_MEFL, 0.90, na.rm = TRUE),
    smoothed_top10 = smoothed_log10_MEFL >= quantile(smoothed_log10_MEFL, 0.90, na.rm = TRUE),
    high_support = bright_total > 5
  )

u_top10_orig <- ggplot(umap_df2, aes(UMAP1, UMAP2, color = bright_top10)) +
  geom_point(alpha = 0.6, size = 1) +
  theme_classic() +
  labs(
    title = "UMAP of ESM2 embeddings",
    subtitle = "Top 10% original brightness highlighted",
    color = "Top 10%"
  )

u_top10_smooth <- ggplot(umap_df2, aes(UMAP1, UMAP2, color = smoothed_top10)) +
  geom_point(alpha = 0.6, size = 1) +
  theme_classic() +
  labs(
    title = "UMAP of ESM2 embeddings",
    subtitle = "Top 10% smoothed brightness highlighted",
    color = "Top 10%"
  )

u_support_bin <- ggplot(umap_df2, aes(UMAP1, UMAP2, color = high_support)) +
  geom_point(alpha = 0.6, size = 1) +
  theme_classic() +
  labs(
    title = "UMAP of ESM2 embeddings",
    subtitle = "Support > 5 highlighted",
    color = "Support > 5"
  )

print(u_top10_orig)
print(u_top10_smooth)
print(u_support_bin)

# ==============================
# Save outputs
# ==============================
ggsave("UMAP_original_log10_MEFL.png", u1, width = 7, height = 5, dpi = 300)
ggsave("UMAP_smoothed_log10_MEFL.png", u2, width = 7, height = 5, dpi = 300)
ggsave("UMAP_support.png", u3, width = 7, height = 5, dpi = 300)

ggsave("UMAP_top10_original.png", u_top10_orig, width = 7, height = 5, dpi = 300)
ggsave("UMAP_top10_smoothed.png", u_top10_smooth, width = 7, height = 5, dpi = 300)
ggsave("UMAP_support_gt5.png", u_support_bin, width = 7, height = 5, dpi = 300)