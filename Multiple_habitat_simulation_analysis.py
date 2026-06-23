# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 17:03:40 2026

@author: hanra
"""

# -*- coding: utf-8 -*-
"""
FINAL CROSS-SCENARIO ANALYSIS FOR POWER-LAW / HABITAT SIMULATION STUDY
========================================================================
This script:
1. Scans all scenario folders
2. Loads each results_clean file
3. Computes method performance across scenarios
4. Identifies best methods by scenario and by true beta
5. Produces publication-quality figures
6. Exports manuscript-ready CSV, XLSX, and LaTeX tables

Author: Your Name
Date: 2026

Tested for folder names such as:
    habitat_visualization_output_R_40_20_n72_12
where:
    sigma_x = 0.40
    sigma_y = 0.20
    n_per_habitat = 72
    n_habitats = 12
"""

import os
import re
import math
import warnings
from pathlib import Path
from typing import Optional, Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

warnings.filterwarnings("ignore")

# ============================================================================
# USER SETTINGS
# ============================================================================

BASE_DIR = Path(r"E:/LMEM_2026/Mevan")   # <-- CHANGE IF NEEDED
OUTPUT_DIR = BASE_DIR / "FINAL_ANALYSIS_MANUSCRIPT"

TOP_N_METHODS_FOR_HEATMAP = 12
TOP_N_METHODS_FOR_SAMPLESIZE_PLOT = 8
MIN_NON_MISSING_PER_METHOD = 5

SAVE_PDF = True
SAVE_PNG = True
DPI = 600

# Matplotlib publication style
plt.rcParams.update({
    "figure.dpi": DPI,
    "savefig.dpi": DPI,
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "axes.linewidth": 0.8,
    "lines.linewidth": 2.0,
    "patch.linewidth": 0.8,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "DejaVu Sans"
})

# ============================================================================
# HELPERS
# ============================================================================

def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def sanitize_filename(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_\-\.]+", "_", s)


def parse_scenario_name(folder_name: str) -> Dict[str, Optional[float]]:
    """
    Parses folder names like:
        habitat_visualization_output_R_40_20_n72_12
        habitat_visualization_output_R_30_30_n36_6
    """
    meta = {
        "scenario": folder_name,
        "sigma_x": np.nan,
        "sigma_y": np.nan,
        "n_per_habitat": np.nan,
        "n_habitats": np.nan,
        "n_total": np.nan,
    }

    # Main expected pattern
    m = re.search(r"_([A-Za-z])_(\d+)_(\d+)_n(\d+)_(\d+)$", folder_name)
    if m:
        sigma_x = float(m.group(2)) / 100.0
        sigma_y = float(m.group(3)) / 100.0
        n_per_habitat = int(m.group(4))
        n_habitats = int(m.group(5))

        meta["sigma_x"] = sigma_x
        meta["sigma_y"] = sigma_y
        meta["n_per_habitat"] = n_per_habitat
        meta["n_habitats"] = n_habitats
        meta["n_total"] = n_per_habitat * n_habitats
        return meta

    # Fallback pattern
    m2 = re.search(r"_(\d+)_(\d+)_n(\d+)_(\d+)$", folder_name)
    if m2:
        sigma_x = float(m2.group(1)) / 100.0
        sigma_y = float(m2.group(2)) / 100.0
        n_per_habitat = int(m2.group(3))
        n_habitats = int(m2.group(4))

        meta["sigma_x"] = sigma_x
        meta["sigma_y"] = sigma_y
        meta["n_per_habitat"] = n_per_habitat
        meta["n_habitats"] = n_habitats
        meta["n_total"] = n_per_habitat * n_habitats
        return meta

    return meta


def find_results_file(folder: Path) -> Optional[Path]:
    """
    Looks for results_clean in common formats.
    """
    candidates = []
    for pattern in [
        "results_clean.csv",
        "results_clean.tsv",
        "results_clean.txt",
        "results_clean.xlsx",
        "results_clean.xls",
        "results_clean.ods",
        "results_clean"
    ]:
        candidates.extend(folder.glob(pattern))

    if not candidates:
        # more permissive fallback
        candidates.extend(folder.glob("results_clean*"))

    return candidates[0] if candidates else None


def load_results_file(path: Path) -> pd.DataFrame:
    """
    Load results_clean from various formats.
    """
    suffix = path.suffix.lower()

    if suffix == ".csv":
        return pd.read_csv(path)
    elif suffix == ".tsv":
        return pd.read_csv(path, sep="\t")
    elif suffix == ".txt":
        # try comma first, then tab
        try:
            return pd.read_csv(path)
        except Exception:
            return pd.read_csv(path, sep="\t")
    elif suffix in [".xlsx", ".xls"]:
        return pd.read_excel(path)
    elif suffix == ".ods":
        return pd.read_excel(path, engine="odf")
    else:
        # last resort
        try:
            return pd.read_csv(path)
        except Exception:
            try:
                return pd.read_excel(path)
            except Exception as e:
                raise ValueError(f"Could not load file: {path}") from e


def sem(x: pd.Series) -> float:
    x = pd.Series(x).dropna()
    if len(x) <= 1:
        return np.nan
    return x.std(ddof=1) / np.sqrt(len(x))


def ci95_mean(x: pd.Series) -> Tuple[float, float]:
    x = pd.Series(x).dropna()
    if len(x) <= 1:
        return (np.nan, np.nan)
    m = x.mean()
    s = x.std(ddof=1)
    h = 1.96 * s / np.sqrt(len(x))
    return (m - h, m + h)


def format_float(x, digits=3):
    if pd.isna(x):
        return ""
    return f"{x:.{digits}f}"


def latex_escape(s: str) -> str:
    if s is None:
        return ""
    s = str(s)
    rep = {
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
        "\\": r"\textbackslash{}",
    }
    for k, v in rep.items():
        s = s.replace(k, v)
    return s


def save_fig(fig: plt.Figure, outbase: Path):
    if SAVE_PNG:
        fig.savefig(str(outbase.with_suffix(".png")))
    if SAVE_PDF:
        fig.savefig(str(outbase.with_suffix(".pdf")))
    plt.close(fig)


def method_pretty_name(m: str) -> str:
    replacements = {
        "ols": "OLS",
        "sma": "SMA",
        "majoraxis": "Major Axis",
        "rma": "RMA",
        "odr": "ODR",
        "deming_std": "Deming",
        "deming_wtd": "Weighted Deming",
        "deming_mcr": "Deming (mcr)",
        "pbablok": "Passing-Bablok",
        "theilsen": "Theil-Sen",
        "siegel": "Siegel",
        "simex": "SIMEX",
        "lmer": "LMM RI",
        "lmer_slopes": "LMM RS",
        "nlme": "NLS",
        "brms_std": "Bayes Std",
        "brms_me": "Bayes ME",
        "brms_robust": "Bayes Robust",
        "brms_horseshoe": "Bayes HS",
        "xgboost": "XGBoost",
        "xgboost_corrected": "XGBoost Corr",
        "rf": "Random Forest",
        "glmnet": "GLMNET",
        "bdcocolasso": "BDCoCoLasso",
        "nnet": "NNET",
        "deepreg": "DeepRegression",
        "keras": "Keras NN",
        "keras_denoised": "Keras DAE",
        "ensemble": "Ensemble"
    }
    return replacements.get(m, m)


# ============================================================================
# DATA INGESTION
# ============================================================================

ensure_dir(OUTPUT_DIR)

scenario_folders = sorted([
    p for p in BASE_DIR.glob("habitat_visualization_output_*")
    if p.is_dir()
])

if len(scenario_folders) == 0:
    raise FileNotFoundError(f"No scenario folders found in {BASE_DIR}")

print(f"Found {len(scenario_folders)} scenario folders.")

all_estimates_long = []
scenario_inventory = []
skipped = []

for folder in scenario_folders:
    scenario_name = folder.name
    meta = parse_scenario_name(scenario_name)
    results_file = find_results_file(folder)

    if results_file is None:
        skipped.append((scenario_name, "No results_clean file found"))
        continue

    try:
        df = load_results_file(results_file)
    except Exception as e:
        skipped.append((scenario_name, f"Failed to read file: {e}"))
        continue

    if "beta_true" not in df.columns:
        skipped.append((scenario_name, "beta_true column not found"))
        continue

    # keep only numeric columns where possible
    for c in df.columns:
        try:
            df[c] = pd.to_numeric(df[c], errors="ignore")
        except Exception:
            pass

    beta_cols = [c for c in df.columns if c.startswith("beta_") and c != "beta_true"]

    if len(beta_cols) == 0:
        skipped.append((scenario_name, "No beta_* method columns found"))
        continue

    cover_cols = [c for c in df.columns if c.startswith("cover_")]
    se_cols = [c for c in df.columns if c.startswith("se_")]

    scenario_inventory.append({
        "scenario": scenario_name,
        "file": str(results_file),
        "n_rows": len(df),
        "n_beta_methods": len(beta_cols),
        "n_cover_cols": len(cover_cols),
        "n_se_cols": len(se_cols),
        **meta
    })

    id_vars = [c for c in ["sim_id", "beta_true", "n", "sigma_x_me", "sigma_y_me",
                           "n_per_habitat", "n_habitats", "time_total"] if c in df.columns]

    long_df = df.melt(
        id_vars=id_vars if len(id_vars) > 0 else ["beta_true"],
        value_vars=beta_cols,
        var_name="method_col",
        value_name="estimate"
    )
    long_df["method"] = long_df["method_col"].str.replace("^beta_", "", regex=True)
    long_df["scenario"] = scenario_name
    long_df["sigma_x"] = meta["sigma_x"]
    long_df["sigma_y"] = meta["sigma_y"]
    long_df["n_per_habitat_scen"] = meta["n_per_habitat"]
    long_df["n_habitats_scen"] = meta["n_habitats"]
    long_df["n_total_scen"] = meta["n_total"]

    # attach coverage if available
    for cc in cover_cols:
        mname = cc.replace("cover_", "")
        mask = long_df["method"] == mname
        if mask.any():
            rep = df[cc].reset_index(drop=True)
            long_df.loc[mask, "coverage"] = np.tile(rep.values, mask.sum() // len(rep) + 1)[:mask.sum()]

    all_estimates_long.append(long_df)

if len(all_estimates_long) == 0:
    raise RuntimeError("No valid scenario results were loaded.")

est = pd.concat(all_estimates_long, ignore_index=True)
inventory_df = pd.DataFrame(scenario_inventory)
skipped_df = pd.DataFrame(skipped, columns=["scenario", "reason"]) if skipped else pd.DataFrame(columns=["scenario", "reason"])

# ============================================================================
# METRICS
# ============================================================================

# Clean estimates
est["estimate"] = pd.to_numeric(est["estimate"], errors="coerce")
est["beta_true"] = pd.to_numeric(est["beta_true"], errors="coerce")

# Remove pathological values
est.loc[np.abs(est["estimate"]) > 1e6, "estimate"] = np.nan

est["error"] = est["estimate"] - est["beta_true"]
est["abs_error"] = np.abs(est["error"])
est["sq_error"] = est["error"] ** 2
est["success"] = (~est["estimate"].isna()).astype(int)

# Scenario x beta_true x method
metrics_sbm = (
    est.groupby(["scenario", "sigma_x", "sigma_y", "n_per_habitat_scen", "n_habitats_scen",
                 "n_total_scen", "beta_true", "method"], dropna=False)
    .agg(
        n_total_obs=("estimate", "size"),
        n_nonmissing=("estimate", lambda x: x.notna().sum()),
        success_rate=("success", "mean"),
        mean_estimate=("estimate", "mean"),
        median_estimate=("estimate", "median"),
        bias=("error", "mean"),
        abs_bias=("error", lambda x: np.abs(np.nanmean(x))),
        mae=("abs_error", "mean"),
        rmse=("sq_error", lambda x: np.sqrt(np.nanmean(x))),
        sd_estimate=("estimate", "std"),
        q25=("estimate", lambda x: np.nanquantile(x.dropna(), 0.25) if x.notna().sum() > 0 else np.nan),
        q75=("estimate", lambda x: np.nanquantile(x.dropna(), 0.75) if x.notna().sum() > 0 else np.nan),
        coverage=("coverage", "mean") if "coverage" in est.columns else ("success", lambda x: np.nan),
    )
    .reset_index()
)

metrics_sbm = metrics_sbm[metrics_sbm["n_nonmissing"] >= MIN_NON_MISSING_PER_METHOD].copy()

# Scenario x method, averaging over beta_true
metrics_sm = (
    metrics_sbm.groupby(["scenario", "sigma_x", "sigma_y", "n_per_habitat_scen",
                         "n_habitats_scen", "n_total_scen", "method"], dropna=False)
    .agg(
        n_beta_levels=("beta_true", "nunique"),
        avg_bias=("bias", "mean"),
        avg_abs_bias=("abs_bias", "mean"),
        avg_mae=("mae", "mean"),
        avg_rmse=("rmse", "mean"),
        avg_coverage=("coverage", "mean"),
        avg_success_rate=("success_rate", "mean"),
        mean_estimate=("mean_estimate", "mean")
    )
    .reset_index()
)

# Method summary over all scenarios
metrics_m = (
    metrics_sm.groupby("method", dropna=False)
    .agg(
        n_scenarios=("scenario", "nunique"),
        mean_rmse=("avg_rmse", "mean"),
        sd_rmse=("avg_rmse", "std"),
        median_rmse=("avg_rmse", "median"),
        mean_mae=("avg_mae", "mean"),
        mean_bias=("avg_bias", "mean"),
        mean_abs_bias=("avg_abs_bias", "mean"),
        mean_coverage=("avg_coverage", "mean"),
        mean_success_rate=("avg_success_rate", "mean")
    )
    .reset_index()
    .sort_values(["mean_rmse", "mean_abs_bias", "mean_mae"], ascending=[True, True, True])
)

metrics_m["rmse_ci_low"], metrics_m["rmse_ci_high"] = zip(*metrics_sm.groupby("method")["avg_rmse"].apply(ci95_mean).reindex(metrics_m["method"]).tolist())
metrics_m["pretty_method"] = metrics_m["method"].map(method_pretty_name)

# Best by scenario
best_by_scenario = (
    metrics_sm.sort_values(
        ["scenario", "avg_rmse", "avg_abs_bias", "avg_mae", "avg_success_rate"],
        ascending=[True, True, True, True, False]
    )
    .groupby("scenario", as_index=False)
    .first()
)

# Best by scenario and beta_true
best_by_scenario_beta = (
    metrics_sbm.sort_values(
        ["scenario", "beta_true", "rmse", "abs_bias", "mae", "success_rate"],
        ascending=[True, True, True, True, True, False]
    )
    .groupby(["scenario", "beta_true"], as_index=False)
    .first()
)

# Win counts
win_counts = (
    best_by_scenario["method"]
    .value_counts()
    .rename_axis("method")
    .reset_index(name="n_wins")
)

win_counts["pretty_method"] = win_counts["method"].map(method_pretty_name)

# Best by sample size and beta_true
best_by_sample_beta = (
    metrics_sbm.sort_values(
        ["n_total_scen", "beta_true", "rmse", "abs_bias", "mae"],
        ascending=[True, True, True, True, True]
    )
    .groupby(["n_total_scen", "beta_true"], as_index=False)
    .first()
)

# Best by noise pair and beta_true
best_by_noise_beta = (
    metrics_sbm.sort_values(
        ["sigma_x", "sigma_y", "beta_true", "rmse", "abs_bias", "mae"],
        ascending=[True, True, True, True, True, True]
    )
    .groupby(["sigma_x", "sigma_y", "beta_true"], as_index=False)
    .first()
)

# ============================================================================
# SAVE RAW ANALYTICS TABLES
# ============================================================================

inventory_df.to_csv(OUTPUT_DIR / "scenario_inventory.csv", index=False)
skipped_df.to_csv(OUTPUT_DIR / "skipped_scenarios.csv", index=False)
est.to_csv(OUTPUT_DIR / "all_estimates_long.csv", index=False)
metrics_sbm.to_csv(OUTPUT_DIR / "metrics_by_scenario_beta_method.csv", index=False)
metrics_sm.to_csv(OUTPUT_DIR / "metrics_by_scenario_method.csv", index=False)
metrics_m.to_csv(OUTPUT_DIR / "overall_method_ranking.csv", index=False)
best_by_scenario.to_csv(OUTPUT_DIR / "best_method_by_scenario.csv", index=False)
best_by_scenario_beta.to_csv(OUTPUT_DIR / "best_method_by_scenario_and_beta.csv", index=False)
best_by_sample_beta.to_csv(OUTPUT_DIR / "best_method_by_sample_size_and_beta.csv", index=False)
best_by_noise_beta.to_csv(OUTPUT_DIR / "best_method_by_noise_and_beta.csv", index=False)
win_counts.to_csv(OUTPUT_DIR / "method_win_counts.csv", index=False)

# Excel workbook
with pd.ExcelWriter(OUTPUT_DIR / "final_manuscript_tables.xlsx", engine="openpyxl") as writer:
    inventory_df.to_excel(writer, sheet_name="inventory", index=False)
    skipped_df.to_excel(writer, sheet_name="skipped", index=False)
    metrics_m.to_excel(writer, sheet_name="overall_ranking", index=False)
    best_by_scenario.to_excel(writer, sheet_name="best_by_scenario", index=False)
    best_by_scenario_beta.to_excel(writer, sheet_name="best_by_scen_beta", index=False)
    best_by_sample_beta.to_excel(writer, sheet_name="best_by_sample_beta", index=False)
    best_by_noise_beta.to_excel(writer, sheet_name="best_by_noise_beta", index=False)
    metrics_sm.to_excel(writer, sheet_name="scenario_method_metrics", index=False)
    metrics_sbm.to_excel(writer, sheet_name="scenario_beta_method", index=False)
    win_counts.to_excel(writer, sheet_name="win_counts", index=False)

print("Saved core tables.")

# ============================================================================
# MANUSCRIPT-READY TABLES
# ============================================================================

# Table 1: overall ranking
table1 = metrics_m.copy()
table1 = table1[[
    "pretty_method", "n_scenarios", "mean_rmse", "sd_rmse", "median_rmse",
    "mean_mae", "mean_bias", "mean_abs_bias", "mean_coverage", "mean_success_rate"
]].rename(columns={
    "pretty_method": "Method",
    "n_scenarios": "Scenarios",
    "mean_rmse": "Mean RMSE",
    "sd_rmse": "SD RMSE",
    "median_rmse": "Median RMSE",
    "mean_mae": "Mean MAE",
    "mean_bias": "Mean Bias",
    "mean_abs_bias": "Mean |Bias|",
    "mean_coverage": "Coverage",
    "mean_success_rate": "Success Rate"
})

# Table 2: best method by scenario
table2 = best_by_scenario.copy()
table2["Method"] = table2["method"].map(method_pretty_name)
table2 = table2[[
    "scenario", "sigma_x", "sigma_y", "n_per_habitat_scen", "n_habitats_scen",
    "n_total_scen", "Method", "avg_rmse", "avg_mae", "avg_bias", "avg_success_rate", "avg_coverage"
]].rename(columns={
    "scenario": "Scenario",
    "sigma_x": "SigmaX",
    "sigma_y": "SigmaY",
    "n_per_habitat_scen": "n per habitat",
    "n_habitats_scen": "Habitats",
    "n_total_scen": "Total n",
    "avg_rmse": "Avg RMSE",
    "avg_mae": "Avg MAE",
    "avg_bias": "Avg Bias",
    "avg_success_rate": "Success Rate",
    "avg_coverage": "Coverage"
})

# Table 3: win counts
table3 = win_counts.copy()
table3 = table3[["pretty_method", "n_wins"]].rename(columns={
    "pretty_method": "Method",
    "n_wins": "Scenario Wins"
})

table1.to_csv(OUTPUT_DIR / "Table1_overall_ranking.csv", index=False)
table2.to_csv(OUTPUT_DIR / "Table2_best_by_scenario.csv", index=False)
table3.to_csv(OUTPUT_DIR / "Table3_win_counts.csv", index=False)

# LaTeX export
def dataframe_to_latex(df: pd.DataFrame, outpath: Path, caption: str, label: str,
                       float_fmt: Dict[str, int] = None):
    df2 = df.copy()
    float_fmt = float_fmt or {}
    for col in df2.columns:
        if pd.api.types.is_numeric_dtype(df2[col]):
            digits = float_fmt.get(col, 3)
            df2[col] = df2[col].map(lambda x: format_float(x, digits))
        else:
            df2[col] = df2[col].map(latex_escape)

    lines = []
    lines.append(r"\begin{table}[!htbp]")
    lines.append(r"\centering")
    lines.append(r"\caption{" + latex_escape(caption) + r"}")
    lines.append(r"\label{" + latex_escape(label) + r"}")
    lines.append(r"\begin{tabular}{" + "l" * len(df2.columns) + r"}")
    lines.append(r"\hline")
    lines.append(" & ".join([latex_escape(c) for c in df2.columns]) + r" \\")
    lines.append(r"\hline")
    for _, row in df2.iterrows():
        lines.append(" & ".join(map(str, row.values)) + r" \\")
    lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    with open(outpath, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

dataframe_to_latex(
    table1, OUTPUT_DIR / "Table1_overall_ranking.tex",
    caption="Overall method ranking across all scenarios based on average RMSE.",
    label="tab:overall_ranking",
    float_fmt={
        "Mean RMSE": 3, "SD RMSE": 3, "Median RMSE": 3,
        "Mean MAE": 3, "Mean Bias": 3, "Mean |Bias|": 3,
        "Coverage": 3, "Success Rate": 3
    }
)

dataframe_to_latex(
    table2, OUTPUT_DIR / "Table2_best_by_scenario.tex",
    caption="Best-performing method within each scenario, ranked by average RMSE across true beta values.",
    label="tab:best_by_scenario",
    float_fmt={
        "SigmaX": 2, "SigmaY": 2, "Avg RMSE": 3, "Avg MAE": 3,
        "Avg Bias": 3, "Success Rate": 3, "Coverage": 3
    }
)

dataframe_to_latex(
    table3, OUTPUT_DIR / "Table3_win_counts.tex",
    caption="Number of scenarios won by each method.",
    label="tab:win_counts",
    float_fmt={"Scenario Wins": 0}
)

print("Saved manuscript-ready tables.")

# ============================================================================
# FIGURES
# ============================================================================

# Add pretty names
metrics_sm["pretty_method"] = metrics_sm["method"].map(method_pretty_name)
metrics_sbm["pretty_method"] = metrics_sbm["method"].map(method_pretty_name)
best_by_scenario["pretty_method"] = best_by_scenario["method"].map(method_pretty_name)
best_by_scenario_beta["pretty_method"] = best_by_scenario_beta["method"].map(method_pretty_name)

# Order methods by overall ranking
method_order = metrics_m["method"].tolist()
pretty_method_order = [method_pretty_name(m) for m in method_order]

# ----------------------------------------------------------------------------
# FIGURE 1: Overall ranking by mean RMSE with 95% CI
# ----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, max(5, 0.35 * len(metrics_m))))
plot_df = metrics_m.copy().sort_values("mean_rmse", ascending=True)

y = np.arange(len(plot_df))
ax.errorbar(
    plot_df["mean_rmse"], y,
    xerr=[
        plot_df["mean_rmse"] - plot_df["rmse_ci_low"],
        plot_df["rmse_ci_high"] - plot_df["mean_rmse"]
    ],
    fmt="o", capsize=3
)
ax.set_yticks(y)
ax.set_yticklabels(plot_df["pretty_method"])
ax.set_xlabel("Mean RMSE across scenarios")
ax.set_ylabel("Method")
ax.set_title("Overall method ranking across scenarios")
ax.invert_yaxis()
ax.grid(False)
fig.tight_layout()
save_fig(fig, OUTPUT_DIR / "Figure1_overall_method_ranking")

# ----------------------------------------------------------------------------
# FIGURE 2: Win counts by method
# ----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 5))
wc = win_counts.copy()
wc = wc.sort_values("n_wins", ascending=False)
ax.bar(wc["pretty_method"], wc["n_wins"])
ax.set_ylabel("Number of scenario wins")
ax.set_xlabel("Method")
ax.set_title("Scenario wins by method")
ax.tick_params(axis="x", rotation=45)
ax.grid(False)
fig.tight_layout()
save_fig(fig, OUTPUT_DIR / "Figure2_method_win_counts")

# ----------------------------------------------------------------------------
# FIGURE 3: Heatmap of average RMSE by scenario and method
# ----------------------------------------------------------------------------
top_methods = metrics_m.head(TOP_N_METHODS_FOR_HEATMAP)["method"].tolist()
heat_df = metrics_sm[metrics_sm["method"].isin(top_methods)].copy()

# scenario ordering
scenario_order_df = (
    inventory_df[["scenario", "sigma_x", "sigma_y", "n_total"]]
    .drop_duplicates()
    .sort_values(["n_total", "sigma_x", "sigma_y", "scenario"])
)
scenario_order = scenario_order_df["scenario"].tolist()

heat_pivot = heat_df.pivot_table(
    index="scenario", columns="method", values="avg_rmse", aggfunc="mean"
).reindex(index=scenario_order, columns=top_methods)

fig, ax = plt.subplots(figsize=(10, max(6, 0.35 * len(heat_pivot))))
im = ax.imshow(heat_pivot.values, aspect="auto")
ax.set_xticks(np.arange(len(top_methods)))
ax.set_xticklabels([method_pretty_name(m) for m in top_methods], rotation=45, ha="right")
ax.set_yticks(np.arange(len(heat_pivot.index)))
ax.set_yticklabels(heat_pivot.index)
ax.set_title("Average RMSE by scenario and method")
ax.set_xlabel("Method")
ax.set_ylabel("Scenario")
cbar = fig.colorbar(im, ax=ax)
cbar.set_label("Average RMSE")
fig.tight_layout()
save_fig(fig, OUTPUT_DIR / "Figure3_heatmap_avg_rmse")

# ----------------------------------------------------------------------------
# FIGURE 4: Mean RMSE versus sample size for top methods
# ----------------------------------------------------------------------------
top_methods_sample = metrics_m.head(TOP_N_METHODS_FOR_SAMPLESIZE_PLOT)["method"].tolist()
sample_df = (
    metrics_sm[metrics_sm["method"].isin(top_methods_sample)]
    .groupby(["n_total_scen", "method"], as_index=False)
    .agg(mean_rmse=("avg_rmse", "mean"))
)

fig, ax = plt.subplots(figsize=(8, 5.5))
for m in top_methods_sample:
    sub = sample_df[sample_df["method"] == m].sort_values("n_total_scen")
    ax.plot(sub["n_total_scen"], sub["mean_rmse"], marker="o", label=method_pretty_name(m))
ax.set_xlabel("Total sample size")
ax.set_ylabel("Mean RMSE")
ax.set_title("Performance versus sample size")
ax.grid(False)
ax.legend(frameon=False, ncol=2)
fig.tight_layout()
save_fig(fig, OUTPUT_DIR / "Figure4_performance_vs_sample_size")

# ----------------------------------------------------------------------------
# FIGURE 5: Bias versus RMSE scatter
# ----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 6))
sc = ax.scatter(metrics_m["mean_abs_bias"], metrics_m["mean_rmse"])
for _, row in metrics_m.iterrows():
    ax.text(row["mean_abs_bias"], row["mean_rmse"], method_pretty_name(row["method"]), fontsize=8)
ax.set_xlabel("Mean absolute bias")
ax.set_ylabel("Mean RMSE")
ax.set_title("Bias-accuracy tradeoff across methods")
ax.grid(False)
fig.tight_layout()
save_fig(fig, OUTPUT_DIR / "Figure5_bias_rmse_tradeoff")

# ----------------------------------------------------------------------------
# FIGURE 6: Winning method landscape by sample size and true beta
# ----------------------------------------------------------------------------
land = best_by_sample_beta.copy()
land["pretty_method"] = land["method"].map(method_pretty_name)

method_levels = metrics_m["method"].tolist()
method_to_int = {m: i for i, m in enumerate(method_levels)}
int_to_label = {i: method_pretty_name(m) for i, m in enumerate(method_levels)}

land["method_id"] = land["method"].map(method_to_int)

pivot_land = land.pivot_table(
    index="beta_true", columns="n_total_scen", values="method_id", aggfunc="first"
).sort_index()

fig, ax = plt.subplots(figsize=(8, 4.5))
cmap = plt.cm.get_cmap("tab20", len(method_levels))
im = ax.imshow(pivot_land.values, aspect="auto", cmap=cmap, interpolation="nearest")
ax.set_xticks(np.arange(len(pivot_land.columns)))
ax.set_xticklabels([str(int(x)) if not pd.isna(x) else "" for x in pivot_land.columns], rotation=0)
ax.set_yticks(np.arange(len(pivot_land.index)))
ax.set_yticklabels([format_float(x, 2) for x in pivot_land.index])
ax.set_xlabel("Total sample size")
ax.set_ylabel(r"True $\beta$")
ax.set_title("Best-performing method by sample size and true beta")
ax.grid(False)

legend_handles = [Patch(facecolor=cmap(i), edgecolor="none", label=int_to_label[i]) for i in sorted(set(land["method_id"].dropna().astype(int)))]
ax.legend(handles=legend_handles, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
fig.tight_layout()
save_fig(fig, OUTPUT_DIR / "Figure6_winning_method_landscape")

# ----------------------------------------------------------------------------
# FIGURE 7: Winning method landscape by noise pair and true beta
# ----------------------------------------------------------------------------
noise = best_by_noise_beta.copy()
noise["noise_label"] = noise.apply(lambda r: f"({r['sigma_x']:.2f}, {r['sigma_y']:.2f})", axis=1)
noise["method_id"] = noise["method"].map(method_to_int)

pivot_noise = noise.pivot_table(
    index="beta_true", columns="noise_label", values="method_id", aggfunc="first"
).sort_index(axis=1)

fig, ax = plt.subplots(figsize=(9, 4.5))
im = ax.imshow(pivot_noise.values, aspect="auto", cmap=cmap, interpolation="nearest")
ax.set_xticks(np.arange(len(pivot_noise.columns)))
ax.set_xticklabels(pivot_noise.columns, rotation=45, ha="right")
ax.set_yticks(np.arange(len(pivot_noise.index)))
ax.set_yticklabels([format_float(x, 2) for x in pivot_noise.index])
ax.set_xlabel(r"Noise pair $(\sigma_x,\sigma_y)$")
ax.set_ylabel(r"True $\beta$")
ax.set_title("Best-performing method by noise level and true beta")
ax.grid(False)
legend_handles = [Patch(facecolor=cmap(i), edgecolor="none", label=int_to_label[i]) for i in sorted(set(noise["method_id"].dropna().astype(int)))]
ax.legend(handles=legend_handles, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
fig.tight_layout()
save_fig(fig, OUTPUT_DIR / "Figure7_winning_method_by_noise")

print("Saved figures.")

# ============================================================================
# SUMMARY TXT
# ============================================================================

with open(OUTPUT_DIR / "analysis_summary.txt", "w", encoding="utf-8") as f:
    f.write("FINAL CROSS-SCENARIO ANALYSIS SUMMARY\n")
    f.write("===================================\n\n")
    f.write(f"Base directory: {BASE_DIR}\n")
    f.write(f"Scenario folders found: {len(scenario_folders)}\n")
    f.write(f"Scenario folders loaded: {inventory_df.shape[0]}\n")
    f.write(f"Scenario folders skipped: {skipped_df.shape[0]}\n\n")

    f.write("Top 10 methods by overall mean RMSE:\n")
    top10 = metrics_m.head(10)
    for i, (_, row) in enumerate(top10.iterrows(), start=1):
        f.write(
            f"{i:2d}. {method_pretty_name(row['method']):<20s}  "
            f"Mean RMSE={row['mean_rmse']:.4f}, "
            f"Mean |Bias|={row['mean_abs_bias']:.4f}, "
            f"Mean MAE={row['mean_mae']:.4f}, "
            f"Scenarios={int(row['n_scenarios'])}\n"
        )

    f.write("\nMethods winning the most scenarios:\n")
    for _, row in win_counts.iterrows():
        f.write(f"  {method_pretty_name(row['method'])}: {int(row['n_wins'])}\n")

print("\nDONE.")
print(f"All outputs saved in:\n{OUTPUT_DIR}")