#!/usr/bin/env python3

import FlowCal
import numpy as np
import pandas as pd
import os

# ============================================================
# SETTINGS
# ============================================================
bead_file = "Samples_Rainbow bead 50000.fcs"

sample_files = [
    "Samples_Brighter blue library_sort.fcs",
    "Samples_Med-Brighter blue library_sort.fcs",
    "Samples_Medium blue library_sort.fcs",
]

blue_channel = "405-1  431_28-A"

# Manufacturer-provided MEFL values for the 7 usable bead peaks
mefl_values = np.array([792, 2079, 6588, 16471, 47497, 137049, 271647])

# Sort gate boundaries on the ORIGINAL sort scale
# include lower bound, exclude upper bound
bins = pd.DataFrame({
    "bin": ["P12", "P13", "P14", "P15", "P16", "P17", "P18", "P19", "P20", "P21", "P22"],
    "min": [0, 41, 67, 108, 171, 279, 439, 726, 1288, 2126, 3294],
    "max": [41, 67, 108, 171, 279, 439, 723, 1282, 2123, 3286, 5260],
})

print("RUNNING UPDATED SCRIPT WITH RAW-BIN / MEFL-MEDIAN LOGIC")
print("SCRIPT PATH:", os.path.abspath(__file__))

# ============================================================
# 1. Build transform from bead file
# ============================================================
print("\nLoading beads and building MEFL transform...")
beads = FlowCal.io.FCSData(bead_file)
beads = FlowCal.transform.to_rfi(beads)

out = FlowCal.mef.get_transform_fxn(
    beads,
    mef_values=mefl_values,
    mef_channels=blue_channel,
    plot=True,
    full_output=True
)

to_mefl = out.transform_fxn

print("\nRaw representative bead fluorescence values (RFI):")
for i, val in enumerate(out.statistic["values"][0], start=1):
    print(f"Peak {i}: {val}")

print("\nSelected RFI values used for fitting:")
print(out.selection["rfi"][0])

print("\nSelected MEFL values used for fitting:")
print(out.selection["mef"][0])

# ============================================================
# 2. Helper: load raw sort-scale values + MEFL values
# ============================================================
def load_sample_raw_and_mefl(sample_file, channel, to_mefl):
    print(f"\nProcessing sample: {sample_file}")

    # Raw FCS data for sort-scale bin assignment
    sample_raw = FlowCal.io.FCSData(sample_file)

    # RFI version for MEFL transformation
    sample_rfi = FlowCal.transform.to_rfi(sample_raw)

    # Apply MEFL transform to the blue channel
    sample_mefl = to_mefl(sample_rfi, channels=channel)

    raw_idx = list(sample_raw.channels).index(channel)
    mefl_idx = list(sample_mefl.channels).index(channel)

    raw_vals = np.asarray(sample_raw[:, raw_idx]).flatten()
    mefl_vals = np.asarray(sample_mefl[:, mefl_idx]).flatten()

    n_total = len(raw_vals)

    # Remove NaN in either vector
    valid_mask = (~np.isnan(raw_vals)) & (~np.isnan(mefl_vals))
    raw_vals = raw_vals[valid_mask]
    mefl_vals = mefl_vals[valid_mask]

    n_after_nan = len(raw_vals)

    # Exclude negative MEFL values
    mefl_nonneg_mask = mefl_vals >= 0
    raw_vals = raw_vals[mefl_nonneg_mask]
    mefl_vals = mefl_vals[mefl_nonneg_mask]

    n_after_mefl_filter = len(raw_vals)

    # Raw negatives are okay to leave here; they simply won't fall into bins
    n_raw_negative = np.sum(raw_vals < 0)

    print(f"  Total events: {n_total}")
    print(f"  Removed due to NaN: {n_total - n_after_nan}")
    print(f"  Removed due to negative MEFL: {n_after_nan - n_after_mefl_filter}")
    print(f"  Events kept: {n_after_mefl_filter}")
    print(f"  Raw values below first bin (<0): {n_raw_negative}")

    if len(mefl_vals) > 0:
        print(f"  MEFL min kept: {np.min(mefl_vals)}")
        print(f"  MEFL median kept: {np.median(mefl_vals)}")
        print(f"  MEFL max kept: {np.max(mefl_vals)}")

    return raw_vals, mefl_vals

# ============================================================
# 3. Helper: calculate per-bin MEFL stats using RAW bin assignment
# ============================================================
def calc_bin_stats_from_raw_bins(raw_vals, mefl_vals, bins_df):
    rows = []

    for _, row in bins_df.iterrows():
        lower = row["min"]
        upper = row["max"]

        in_bin = (raw_vals >= lower) & (raw_vals < upper)
        mefl_in_bin = mefl_vals[in_bin]

        rows.append({
            "bin": row["bin"],
            "min": lower,
            "max": upper,
            "n_events": len(mefl_in_bin),
            "median_brightness": np.median(mefl_in_bin) if len(mefl_in_bin) else np.nan,
            "mean_brightness": np.mean(mefl_in_bin) if len(mefl_in_bin) else np.nan,
            "min_observed": np.min(mefl_in_bin) if len(mefl_in_bin) else np.nan,
            "max_observed": np.max(mefl_in_bin) if len(mefl_in_bin) else np.nan,
        })

    return pd.DataFrame(rows)

# ============================================================
# 4. Process each individual sample
# ============================================================
sample_data = {}
all_results = []

for sample_file in sample_files:
    raw_vals, mefl_vals = load_sample_raw_and_mefl(sample_file, blue_channel, to_mefl)
    sample_data[sample_file] = {
        "raw": raw_vals,
        "mefl": mefl_vals
    }

    sample_base = sample_file.replace(".fcs", "").replace(" ", "_")

    # Save event-level table with both raw and MEFL values
    event_df = pd.DataFrame({
        f"{blue_channel}_raw": raw_vals,
        f"{blue_channel}_MEFL": mefl_vals
    })
    event_outfile = f"{sample_base}_raw_and_mefl_events_nonnegativeMEFL.csv"
    event_df.to_csv(event_outfile, index=False)
    print(f"  Wrote event-level values: {event_outfile}")

    # Save per-bin stats: raw-defined bins, MEFL summaries
    bin_stats = calc_bin_stats_from_raw_bins(raw_vals, mefl_vals, bins)
    bin_stats.insert(0, "sample_file", sample_file)

    stats_outfile = f"{sample_base}_bin_stats_mefl_from_rawbins.csv"
    bin_stats.to_csv(stats_outfile, index=False)
    print(f"  Wrote bin stats: {stats_outfile}")

    all_results.append(bin_stats)

# ============================================================
# 5. Combine Brighter + Medium, since those were mixed in NGS
# ============================================================
combined_name = "Samples_Brighter+Medium_combined"

combined_raw = np.concatenate([
    sample_data["Samples_Brighter blue library_sort.fcs"]["raw"],
    sample_data["Samples_Medium blue library_sort.fcs"]["raw"]
])

combined_mefl = np.concatenate([
    sample_data["Samples_Brighter blue library_sort.fcs"]["mefl"],
    sample_data["Samples_Medium blue library_sort.fcs"]["mefl"]
])

print(f"\nProcessing combined sample: {combined_name}")
print(f"  Combined events kept: {len(combined_mefl)}")

combined_event_df = pd.DataFrame({
    f"{blue_channel}_raw": combined_raw,
    f"{blue_channel}_MEFL": combined_mefl
})
combined_event_outfile = "Samples_Brighter_plus_Medium_combined_raw_and_mefl_events_nonnegativeMEFL.csv"
combined_event_df.to_csv(combined_event_outfile, index=False)
print(f"  Wrote combined event-level values: {combined_event_outfile}")

combined_bin_stats = calc_bin_stats_from_raw_bins(combined_raw, combined_mefl, bins)
combined_bin_stats.insert(0, "sample_file", combined_name)

combined_stats_outfile = "Samples_Brighter_plus_Medium_combined_bin_stats_mefl_from_rawbins.csv"
combined_bin_stats.to_csv(combined_stats_outfile, index=False)
print(f"  Wrote combined bin stats: {combined_stats_outfile}")

all_results.append(combined_bin_stats)

# ============================================================
# 6. Save all summaries together
# ============================================================
combined_summary = pd.concat(all_results, ignore_index=True)
combined_summary.to_csv("all_samples_bin_stats_mefl_from_rawbins.csv", index=False)
print("\nWrote combined summary: all_samples_bin_stats_mefl_from_rawbins.csv")