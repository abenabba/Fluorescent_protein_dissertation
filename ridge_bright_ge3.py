import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score, mean_squared_error
from scipy.stats import pearsonr

# --------------------------------------------------
# Paths
# --------------------------------------------------
embeddings_file = "../esm2_embeddings_len150.npy"
meta_file = "esm2_embeddings_len150_metadata.csv"
target_file = "bright_ge3_barcode_clusters_len150.csv"

# --------------------------------------------------
# Load data
# --------------------------------------------------
X = np.load(embeddings_file)
meta = pd.read_csv(meta_file)
targets = pd.read_csv(target_file)

print("Embeddings shape:", X.shape)
print("Metadata shape:", meta.shape)
print("Targets shape:", targets.shape)

# --------------------------------------------------
# Basic checks
# --------------------------------------------------
required_meta = ["aa_seq"]
required_targets = ["aa_seq", "log10_MEFL", "n_barcode_clusters_total", "aa_length"]

for col in required_meta:
    if col not in meta.columns:
        raise ValueError(f"Missing required metadata column: {col}")

for col in required_targets:
    if col not in targets.columns:
        raise ValueError(f"Missing required target column: {col}")

# --------------------------------------------------
# Align embeddings to filtered target set
# --------------------------------------------------
if "embedding_row" in meta.columns:
    merged = meta.merge(
        targets[["aa_seq", "log10_MEFL", "n_barcode_clusters_total", "aa_length"]],
        on="aa_seq",
        how="inner"
    )
    X_aligned = X[merged["embedding_row"].values]
else:
    meta = meta.copy()
    meta["row_id"] = np.arange(len(meta))
    merged = meta.merge(
        targets[["aa_seq", "log10_MEFL", "n_barcode_clusters_total", "aa_length"]],
        on="aa_seq",
        how="inner"
    )
    X_aligned = X[merged["row_id"].values]

print("Merged shape:", merged.shape)
print("Aligned embedding matrix shape:", X_aligned.shape)

# --------------------------------------------------
# Target
# --------------------------------------------------
y = merged["log10_MEFL"].values

mask = np.isfinite(y)
X_use = X_aligned[mask]
y_use = y[mask]

print("Usable rows:", len(y_use))
print("Target summary:")
print(pd.Series(y_use).describe())

# --------------------------------------------------
# Train/test split
# --------------------------------------------------
X_train, X_test, y_train, y_test = train_test_split(
    X_use, y_use, test_size=0.2, random_state=42
)

print("Train size:", len(y_train))
print("Test size :", len(y_test))

# --------------------------------------------------
# Ridge regression
# --------------------------------------------------
model = Ridge(alpha=1.0)
model.fit(X_train, y_train)
y_pred = model.predict(X_test)

# --------------------------------------------------
# Metrics
# --------------------------------------------------
r2 = r2_score(y_test, y_pred)
rmse = np.sqrt(mean_squared_error(y_test, y_pred))
pearson, pearson_p = pearsonr(y_test, y_pred)

print("\nResults for bright_ge3 log10_MEFL")
print("R^2       :", r2)
print("RMSE      :", rmse)
print("Pearson r :", pearson)
print("Pearson p :", pearson_p)

# --------------------------------------------------
# Save predictions
# --------------------------------------------------
results_df = pd.DataFrame({
    "y_test": y_test,
    "y_pred": y_pred
})

results_df.to_csv("ridge_bright_ge3_predictions.csv", index=False)

summary_df = pd.DataFrame([{
    "dataset": "bright_ge3_barcode_clusters_len150",
    "target": "log10_MEFL",
    "n_total": len(y_use),
    "n_train": len(y_train),
    "n_test": len(y_test),
    "r2": r2,
    "rmse": rmse,
    "pearson": pearson,
    "pearson_p": pearson_p
}])

summary_df.to_csv("ridge_bright_ge3_summary.csv", index=False)

print("\nSaved:")
print(" - ridge_bright_ge3_predictions.csv")
print(" - ridge_bright_ge3_summary.csv")