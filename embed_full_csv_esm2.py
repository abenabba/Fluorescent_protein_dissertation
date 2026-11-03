import numpy as np
import pandas as pd
import torch
from transformers import AutoTokenizer, EsmModel

INPUT_CSV = "brighter_flowseq_inferred_brightness_updated.csv"
SEQ_COL = "aa_seq"
TARGET_COL = "median_MEFL"
MODEL_NAME = "facebook/esm2_t30_150M_UR50D"

MIN_LEN = 150
BATCH_SIZE = 8

OUT_EMBED_NPY = "esm2_embeddings_len150.npy"
OUT_META_CSV = "esm2_embeddings_len150_metadata.csv"

device = "cuda" if torch.cuda.is_available() else "cpu"
print("Using device:", device)

# -----------------------------
# Load and clean data
# -----------------------------
df = pd.read_csv(INPUT_CSV)
df = df[[SEQ_COL, TARGET_COL]].dropna().copy()

df[SEQ_COL] = df[SEQ_COL].astype(str).str.strip().str.upper()
df[SEQ_COL] = df[SEQ_COL].str.replace(" ", "", regex=False)

valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
df = df[df[SEQ_COL].apply(lambda s: len(s) > 0 and set(s).issubset(valid_aas))].copy()
df["aa_length"] = df[SEQ_COL].str.len()
df = df[df["aa_length"] >= MIN_LEN].reset_index(drop=True)

print("Sequences retained:", len(df))

# -----------------------------
# Load model
# -----------------------------
tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
model = EsmModel.from_pretrained(MODEL_NAME)
model.to(device)
model.eval()

# -----------------------------
# Mean pooling
# -----------------------------
def mean_pool_hidden(hidden_states, attention_mask):
    mask = attention_mask.unsqueeze(-1).expand(hidden_states.size()).float()
    summed = torch.sum(hidden_states * mask, dim=1)
    counts = torch.clamp(mask.sum(dim=1), min=1e-9)
    return summed / counts

# -----------------------------
# Embed in batches
# -----------------------------
all_embeddings = []

with torch.no_grad():
    for start in range(0, len(df), BATCH_SIZE):
        end = min(start + BATCH_SIZE, len(df))
        seqs = df[SEQ_COL].iloc[start:end].tolist()

        inputs = tokenizer(
            seqs,
            return_tensors="pt",
            padding=True,
            truncation=True
        )
        inputs = {k: v.to(device) for k, v in inputs.items()}

        outputs = model(**inputs)
        pooled = mean_pool_hidden(outputs.last_hidden_state, inputs["attention_mask"])
        pooled = pooled.cpu().numpy()

        all_embeddings.append(pooled)

        print(f"Processed {end}/{len(df)}")

embeddings = np.vstack(all_embeddings)

# -----------------------------
# Save outputs
# -----------------------------
np.save(OUT_EMBED_NPY, embeddings)

meta = df.copy()
meta["embedding_row"] = np.arange(len(meta))
meta.to_csv(OUT_META_CSV, index=False)

print("Done.")
print("Embedding matrix shape:", embeddings.shape)
print("Saved:", OUT_EMBED_NPY)
print("Saved:", OUT_META_CSV)
