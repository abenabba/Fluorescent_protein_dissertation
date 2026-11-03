library(tidyverse)
library(RcppCNPy)

embedding_file <- "/Users/anissab/Desktop/illumina_shuffled_multi-bin_barcode_analysis/af3_analysis/all_brighter_af2_single_repr_meanpooled.npy"
metadata_file <- "/Users/anissab/Desktop/illumina_shuffled_multi-bin_barcode_analysis/af3_analysis/all_brighter_af2_single_repr_metadata.csv"

X <- RcppCNPy::npyLoad(embedding_file)
meta <- read_csv(metadata_file, show_col_types = FALSE)

dim(X)
nrow(meta)
names(meta)

# check row norms
norms <- sqrt(rowSums(X^2))

summary(norms)

plot(norms, pch = 16, cex = 0.4,
     main = "AF2 embedding row norms",
     xlab = "row index",
     ylab = "L2 norm")

# check embedding variance by dimension
dim_var <- apply(X, 2, var)

summary(dim_var)

plot(dim_var, type = "h",
     main = "Variance per AF2 embedding dimension",
     xlab = "embedding dimension",
     ylab = "variance")

# check exact duplicate vectors
dup_n <- sum(duplicated(as.data.frame(round(X, 6))))
cat("Approx duplicate rows:", dup_n, "\n")

# check PCA variance
pca <- prcomp(scale(X), center = FALSE, scale. = FALSE)

var_explained <- pca$sdev^2 / sum(pca$sdev^2)

plot(var_explained[1:50], type = "b",
     main = "PCA variance explained",
     xlab = "PC",
     ylab = "Proportion variance")