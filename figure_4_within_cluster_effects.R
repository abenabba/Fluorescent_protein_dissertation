# =========================================================
# Figure 4.X
# Within-cluster effects of chromophore-local structural
# features on residual fluorescence brightness
# =========================================================

library(tidyverse)
library(patchwork)

# ---------------------------------------------------------
# Paths
# ---------------------------------------------------------
base_dir <- "~/Desktop/illumina_shuffled_multi-bin_barcode_analysis"
analysis_dir <- file.path(base_dir, "af3_analysis")

input_file <- file.path(
  analysis_dir,
  "bright_af2_chromophore_nonfloor_all_support_k4",
  "chromophore_nonfloor_all_support_clustered_k4.csv"
)

out_dir <- file.path(
  analysis_dir,
  "figure_4_within_cluster_effects"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------
# Settings
# ---------------------------------------------------------
font_family <- "Times New Roman"

structural_features <- c(
  "motif_plddt",
  "local_plddt",
  "motif_span",
  "local_span",
  "n_ca_within_radius",
  "mean_ca_dist_within_radius",
  "n_atoms_within_radius",
  "n_aromatic_ca_within_radius",
  "n_polar_ca_within_radius",
  "n_charged_ca_within_radius",
  "aromatic_n",
  "polar_n",
  "charged_n",
  "gly_pro_n"
)

feature_labels <- c(
  motif_plddt = "Chromophore motif pLDDT",
  local_plddt = "Local pLDDT",
  motif_span = "Chromophore motif span",
  local_span = "Local structural span",
  n_ca_within_radius = "Nearby Cα count",
  mean_ca_dist_within_radius = "Mean nearby Cα distance",
  n_atoms_within_radius = "Nearby atom count",
  n_aromatic_ca_within_radius = "Nearby aromatic residues",
  n_polar_ca_within_radius = "Nearby polar residues",
  n_charged_ca_within_radius = "Nearby charged residues",
  aromatic_n = "Aromatic residues in window",
  polar_n = "Polar residues in window",
  charged_n = "Charged residues in window",
  gly_pro_n = "Gly/Pro residues in window"
)

# ---------------------------------------------------------
# Load clustered dataset
# ---------------------------------------------------------
df <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(
    cluster = factor(cluster, levels = c("1", "2", "3", "4"))
  )

cat("Rows loaded:", nrow(df), "\n")
print(table(df$cluster))

# ---------------------------------------------------------
# Recompute residual brightness
# ---------------------------------------------------------
resid_model <- lm(
  log10_MEFL ~ log10(barcode_cluster_support + 1) + aa_length,
  data = df
)

df <- df %>%
  mutate(
    predicted_from_support_length = predict(resid_model, newdata = df),
    residual_brightness = log10_MEFL - predicted_from_support_length
  )

# ---------------------------------------------------------
# Helper: Cohen's d
# ---------------------------------------------------------
cohens_d_fun <- function(x_top, x_bottom) {
  x_top <- x_top[is.finite(x_top)]
  x_bottom <- x_bottom[is.finite(x_bottom)]

  n1 <- length(x_top)
  n2 <- length(x_bottom)

  if (n1 < 2 || n2 < 2) return(NA_real_)

  v1 <- var(x_top)
  v2 <- var(x_bottom)

  if (!is.finite(v1) || !is.finite(v2)) return(NA_real_)

  pooled_sd <- sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))

  if (!is.finite(pooled_sd) || pooled_sd == 0) return(NA_real_)

  (mean(x_top, na.rm = TRUE) - mean(x_bottom, na.rm = TRUE)) / pooled_sd
}

# ---------------------------------------------------------
# Within-cluster feature effects using residual brightness
# ---------------------------------------------------------
feature_effects_residual <- map_dfr(
  levels(df$cluster),
  function(clust) {

    sub <- df %>%
      filter(cluster == clust)

    top_cut <- quantile(sub$residual_brightness, 0.75, na.rm = TRUE)
    bottom_cut <- quantile(sub$residual_brightness, 0.25, na.rm = TRUE)

    extreme <- sub %>%
      filter(
        residual_brightness >= top_cut |
          residual_brightness <= bottom_cut
      ) %>%
      mutate(
        brightness_group = if_else(
          residual_brightness >= top_cut,
          "top",
          "bottom"
        )
      )

    map_dfr(structural_features, function(feat) {

      x_top <- extreme %>%
        filter(brightness_group == "top") %>%
        pull(all_of(feat))

      x_bottom <- extreme %>%
        filter(brightness_group == "bottom") %>%
        pull(all_of(feat))

      tibble(
        cluster = clust,
        cluster_n = nrow(sub),
        feature = feat,
        n_top = sum(is.finite(x_top)),
        n_bottom = sum(is.finite(x_bottom)),
        mean_top = mean(x_top, na.rm = TRUE),
        mean_bottom = mean(x_bottom, na.rm = TRUE),
        cohens_d = cohens_d_fun(x_top, x_bottom),
        p_value = tryCatch(
          t.test(x_top, x_bottom)$p.value,
          error = function(e) NA_real_
        )
      )
    })
  }
) %>%
  group_by(cluster) %>%
  mutate(
    p_BH_within_cluster = p.adjust(p_value, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    p_BH_global = p.adjust(p_value, method = "BH"),
    feature_label = recode(feature, !!!feature_labels)
  )

write_csv(
  feature_effects_residual,
  file.path(out_dir, "within_cluster_feature_effects_residual.csv")
)

# ---------------------------------------------------------
# Panel A: heatmap of Cohen's d across all clusters/features
# ---------------------------------------------------------
heatmap_df <- feature_effects_residual %>%
  mutate(
    feature_label = factor(
      feature_label,
      levels = rev(unname(feature_labels[structural_features]))
    ),
    sig_within_cluster = p_BH_within_cluster < 0.05
  )

pA <- ggplot(
  heatmap_df,
  aes(x = cluster, y = feature_label, fill = cohens_d)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = ifelse(sig_within_cluster, "*", "")),
    size = 7,
    family = font_family,
    fontface = "bold"
  ) +
  scale_fill_gradient2(
    low = "#3B7EA1",
    mid = "white",
    high = "#D95F02",
    midpoint = 0,
    limits = c(-0.8, 0.8),
    breaks = c(-0.8, -0.4, 0, 0.4, 0.8),
    name = "Cohen's d\n(top - bottom)"
  ) +
  labs(
    title = "Within-cluster chromophore-local feature effects",
    subtitle = "* indicates BH-adjusted significance within that cluster only",
    x = "Structural cluster",
    y = NULL
  ) +
  theme_classic(base_size = 16, base_family = font_family) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    plot.margin = margin(10, 10, 10, 10)
  )

# ---------------------------------------------------------
# Panel B: cluster 2 aromatic-residue boxplot
# ---------------------------------------------------------
cluster2_df <- df %>%
  filter(cluster == "2")

top_cut_c2 <- quantile(cluster2_df$residual_brightness, 0.75, na.rm = TRUE)
bottom_cut_c2 <- quantile(cluster2_df$residual_brightness, 0.25, na.rm = TRUE)

cluster2_extreme <- cluster2_df %>%
  filter(
    residual_brightness >= top_cut_c2 |
      residual_brightness <= bottom_cut_c2
  ) %>%
  mutate(
    brightness_group = if_else(
      residual_brightness >= top_cut_c2,
      "Top 25%",
      "Bottom 25%"
    ),
    brightness_group = factor(
      brightness_group,
      levels = c("Bottom 25%", "Top 25%")
    )
  )

cluster2_stats <- feature_effects_residual %>%
  filter(
    cluster == "2",
    feature == "n_aromatic_ca_within_radius"
  )

subtitle_text <- paste0(
  "n = ", cluster2_stats$n_top, " per group; ",
  "Cohen's d = ", sprintf("%.2f", cluster2_stats$cohens_d), "; ",
  "Welch p = ", sprintf("%.4f", cluster2_stats$p_value), "\n",
  "BH within cluster = ", sprintf("%.3f", cluster2_stats$p_BH_within_cluster), "; ",
  "BH global = ", sprintf("%.3f", cluster2_stats$p_BH_global)
)

pB <- ggplot(
  cluster2_extreme,
  aes(
    x = brightness_group,
    y = n_aromatic_ca_within_radius,
    fill = brightness_group
  )
) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.85,
    color = "black"
  ) +
  geom_jitter(
    width = 0.10,
    height = 0,
    size = 2,
    alpha = 0.65,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Bottom 25%" = "#4C78A8",
      "Top 25%" = "#F58518"
    )
  ) +
  labs(
    title = "Cluster 2: nearby aromatic residues within 10 Å",
    subtitle = subtitle_text,
    x = NULL,
    y = "Nearby aromatic residues"
  ) +
  theme_classic(base_size = 16, base_family = font_family) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )

# ---------------------------------------------------------
# Combine panels
# ---------------------------------------------------------
fig_4x <- pA + pB +
  plot_layout(widths = c(1.6, 1.1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(
      family = font_family,
      size = 20,
      face = "bold"
    )
  )

print(fig_4x)

# ---------------------------------------------------------
# Save figure
# ---------------------------------------------------------
ggsave(
  filename = file.path(out_dir, "Figure_4X_within_cluster_feature_effects.png"),
  plot = fig_4x,
  width = 14,
  height = 8,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(out_dir, "Figure_4X_within_cluster_feature_effects.pdf"),
  plot = fig_4x,
  width = 14,
  height = 8,
  bg = "white"
)