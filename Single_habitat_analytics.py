# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 11:43:41 2026

@author: hanra
"""

# Full single-habitat meta-analysis script
# Fully fixed so beta_lgb / beta_lightgbm / beta_nnet are excluded everywhere.

import re
import warnings
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.tree import DecisionTreeClassifier, plot_tree, export_text
from sklearn.preprocessing import LabelEncoder

warnings.filterwarnings("ignore")

# =============================================================================
# USER SETTINGS
# =============================================================================

BASE_DIR = Path(r"E:\LMEM_2026_Mevan\Single_Habitat")
OUT_DIR = BASE_DIR / f"SingleHabitat_MetaAnalysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TRUE_SLOPES_TARGET = [0.6, 1.0, 1.4]
ERROR_RATIOS_TARGET = [0.5, 1.0, 2.0]
N_TARGET = [36, 72]

# =============================================================================
# PLOTTING STYLE
# =============================================================================

plt.style.use("seaborn-v0_8-whitegrid")
sns.set_context("talk")
plt.rcParams["figure.dpi"] = 300
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["font.size"] = 11
plt.rcParams["axes.titlesize"] = 15
plt.rcParams["axes.labelsize"] = 12
plt.rcParams["legend.fontsize"] = 10

# =============================================================================
# METHOD LABELS
# =============================================================================
# LightGBM and NNET intentionally excluded.

METHOD_LABELS = {
    "beta_ols": "OLS",
    "beta_sma": "SMA",
    "beta_majoraxis": "SMA (Major Axis)",
    "beta_rma": "RMA",
    "beta_odr": "ODR",
    "beta_deming": "Deming",
    "beta_deming_std": "Deming (Std)",
    "beta_deming_wtd": "Deming (Wtd)",
    "beta_deming_mcr": "Deming (MCR)",
    "beta_pbablok": "Passing-Bablok",
    "beta_theilsen": "Theil-Sen",
    "beta_siegel": "Siegel",
    "beta_huber": "Huber",
    "beta_quantile": "Quantile",
    "beta_ridge": "Ridge",
    "beta_lasso": "Lasso",
    "beta_enet": "Elastic Net",
    "beta_adaptive_lasso": "Adaptive Lasso",
    "beta_rf": "Random Forest",
    "beta_gb": "Gradient Boosting",
    "beta_xgb": "XGBoost",
    "beta_xgboost": "XGBoost",
    "beta_xgboost_corrected": "XGBoost Corrected",
    "beta_svm_linear": "SVM Linear",
    "beta_svm_rbf": "SVM RBF",
    "beta_gp": "Gaussian Process",
    "beta_clustered": "Clustered SE",
    "beta_lmer_int": "LMM (Random Intercept)",
    "beta_lmer_slope": "LMM (Random Slope)",
    "beta_lmer": "LMM",
    "beta_lmer_slopes": "LMM (Rand. Slopes)",
    "beta_nlme": "NLME",
    "beta_simex": "SIMEX",
    "beta_glmnet": "GLMNET",
    "beta_keras": "Keras",
    "beta_keras_denoised": "Keras (Denoised)",
    "beta_ensemble": "Ensemble",
    "beta_trimmed": "Trimmed Mean",
    "beta_median": "Median",
    "beta_brms_std": "BRMS Std",
    "beta_brms_me": "BRMS ME",
    "beta_brms_robust": "BRMS Robust",
    "beta_brms_horseshoe": "BRMS Horseshoe",
}

# =============================================================================
# EXCLUDED METHODS
# =============================================================================

EXCLUDED_METHOD_COLUMNS = {
    # Already excluded
    "beta_nnet",
    "beta_lgb",
    "beta_lightgbm",

    # NEW: remove tree-based ML
    "beta_rf",
    "beta_gb",
    "beta_xgb",
    "beta_xgboost",
    "beta_xgboost_corrected",
}

EXCLUDED_METHOD_NAMES = {
    "Lgb",
    "Lightgbm",
    "LightGBM",
    "Nnet",
    "NNET",

    # NEW: remove tree-based ML labels
    "Random Forest",
    "Gradient Boosting",
    "XGBoost",
    "XGBoost Corrected",
}

# =============================================================================
# HELPERS
# =============================================================================

def normalize_numeric(x):
    return pd.to_numeric(x, errors="coerce")


def build_method_name(col):
    if col in EXCLUDED_METHOD_COLUMNS:
        return None
    return METHOD_LABELS.get(col, col.replace("beta_", "").replace("_", " ").title())


def closest_exact_filter(df, col, allowed_values, tol=1e-9):
    keep_mask = np.zeros(len(df), dtype=bool)
    vals = normalize_numeric(df[col]).values
    for av in allowed_values:
        keep_mask |= np.isclose(vals, av, atol=tol, rtol=0)
    return df.loc[keep_mask].copy()


def parse_folder_name(folder_name):
    """
    Expected examples:
      habitat_visualization_single_output_n36_mex20_mey10
      habitat_visualization_single_output_n72_mex20_mey40
    """
    patterns = [
        r"n(\d+)_mex(\d+)_mey(\d+)",
        r"n(\d+)_sx(\d+)_sy(\d+)",
        r"n(\d+)_x(\d+)_y(\d+)",
    ]
    for pat in patterns:
        m = re.search(pat, folder_name)
        if m:
            return int(m.group(1)), float(m.group(2)), float(m.group(3))
    return None


def load_all_single_habitat_data(base_dir):
    all_results = []

    if not base_dir.exists():
        raise FileNotFoundError(f"Base directory not found: {base_dir}")

    folders = [
        f for f in base_dir.iterdir()
        if f.is_dir() and "single" in f.name.lower()
    ]

    if not folders:
        raise ValueError(f"No single-habitat folders found under: {base_dir}")

    print(f"Found {len(folders)} candidate scenario folders")

    for folder in folders:
        parsed = parse_folder_name(folder.name)
        if parsed is None:
            print(f"Skipping folder with unrecognized naming pattern: {folder.name}")
            continue

        n_val, x_error, y_error = parsed

        candidate_files = [
            folder / "all_results.csv",
            folder / "results.csv",
            folder / "all_results.xlsx",
        ]

        result_file = None
        for f in candidate_files:
            if f.exists():
                result_file = f
                break

        if result_file is None:
            print(f"Warning: no recognized results file in {folder.name}")
            continue

        try:
            if result_file.suffix.lower() == ".csv":
                df = pd.read_csv(result_file)
            else:
                df = pd.read_excel(result_file)

            df["n"] = n_val
            df["sigma_x_me"] = x_error
            df["sigma_y_me"] = y_error
            df["error_ratio"] = df["sigma_y_me"] / df["sigma_x_me"]
            df["scenario"] = f"n={n_val}, σx={x_error:g}, σy={y_error:g}"

            all_results.append(df)
            print(f"Loaded: {folder.name} -> {len(df)} rows")
        except Exception as e:
            print(f"Error loading {folder.name}: {e}")

    if not all_results:
        raise ValueError("No result files were successfully loaded.")

    return pd.concat(all_results, ignore_index=True)


def prepare_dataframe(df):
    df = df.copy()

    for c in ["beta_true", "sigma_x_me", "sigma_y_me", "n"]:
        if c in df.columns:
            df[c] = normalize_numeric(df[c])

    if "beta_true" not in df.columns:
        raise ValueError("Column 'beta_true' is required.")

    if "n" not in df.columns:
        raise ValueError("Column 'n' is required for single-habitat analysis.")

    if "sigma_x_me" not in df.columns or "sigma_y_me" not in df.columns:
        raise ValueError("Need sigma_x_me and sigma_y_me columns.")

    df["error_ratio"] = df["sigma_y_me"] / df["sigma_x_me"]
    return df


def find_method_columns(df):
    """
    Find beta_* columns, excluding beta_true and explicitly blocked methods.
    """
    method_cols = []

    for col in df.columns:
        if (
            col.startswith("beta_")
            and col != "beta_true"
            and col not in EXCLUDED_METHOD_COLUMNS
        ):
            method_cols.append(col)

    # Keep known methods in label order first
    ordered = [c for c in METHOD_LABELS if c in method_cols]

    # Then any extra beta_* columns not in METHOD_LABELS but not excluded
    extras = [c for c in method_cols if c not in ordered]

    return ordered + extras


def compute_long_metrics(df, method_cols):
    records = []

    for mcol in method_cols:
        if mcol in EXCLUDED_METHOD_COLUMNS:
            continue

        method_name = build_method_name(mcol)
        if method_name is None or method_name in EXCLUDED_METHOD_NAMES:
            continue

        est = normalize_numeric(df[mcol])
        truth = normalize_numeric(df["beta_true"])
        bias = est - truth
        rmse_point = np.sqrt((est - truth) ** 2)

        tmp = pd.DataFrame({
            "method_col": mcol,
            "method": method_name,
            "beta_true": truth,
            "estimate": est,
            "bias": bias,
            "rmse_point": rmse_point,
            "sigma_x_me": normalize_numeric(df["sigma_x_me"]),
            "sigma_y_me": normalize_numeric(df["sigma_y_me"]),
            "error_ratio": normalize_numeric(df["error_ratio"]),
            "n": normalize_numeric(df["n"]),
            "scenario": df["scenario"] if "scenario" in df.columns else None,
        })

        if "sim_id" in df.columns:
            tmp["sim_id"] = df["sim_id"]
        else:
            tmp["sim_id"] = np.arange(len(df))

        records.append(tmp)

    if not records:
        raise ValueError("No valid method columns left after exclusions.")

    long_df = pd.concat(records, ignore_index=True)
    long_df = long_df.dropna(subset=["estimate", "beta_true"])

    # Final safety filter
    long_df = long_df[~long_df["method_col"].isin(EXCLUDED_METHOD_COLUMNS)].copy()
    long_df = long_df[~long_df["method"].isin(EXCLUDED_METHOD_NAMES)].copy()

    return long_df


def summarize_metrics(long_df):
    summary = (
        long_df
        .groupby(
            ["method", "method_col", "beta_true", "sigma_x_me", "sigma_y_me", "error_ratio", "n"],
            dropna=False,
        )
        .agg(
            mean_bias=("bias", "mean"),
            rmse=("rmse_point", "mean"),
            sd_bias=("bias", "std"),
            sd_rmse=("rmse_point", "std"),
            n_sim=("bias", "size"),
        )
        .reset_index()
    )

    summary = summary[~summary["method_col"].isin(EXCLUDED_METHOD_COLUMNS)].copy()
    summary = summary[~summary["method"].isin(EXCLUDED_METHOD_NAMES)].copy()

    return summary


def save_table(df, name):
    path = OUT_DIR / f"{name}.csv"
    df.to_csv(path, index=False)
    print(f"Saved: {path}")


def save_excel(tables):
    path = OUT_DIR / "summary_tables.xlsx"
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for sheet, table in tables.items():
            table.to_excel(writer, sheet_name=sheet[:31], index=False)
    print(f"Saved: {path}")

# =============================================================================
# PLOTS
# =============================================================================

def plot_bias_rmse_heatmaps_exact_slopes(summary_df):
    plot_df = closest_exact_filter(summary_df, "beta_true", TRUE_SLOPES_TARGET)

    bias_pivot = (
        plot_df.groupby(["method", "beta_true"])["mean_bias"]
        .mean()
        .unstack()
        .reindex(columns=TRUE_SLOPES_TARGET)
    )

    rmse_pivot = (
        plot_df.groupby(["method", "beta_true"])["rmse"]
        .mean()
        .unstack()
        .reindex(columns=TRUE_SLOPES_TARGET)
    )

    bias_pivot = bias_pivot.dropna(how="all", axis=0)
    rmse_pivot = rmse_pivot.dropna(how="all", axis=0)

    method_order = rmse_pivot.mean(axis=1).sort_values().index
    bias_pivot = bias_pivot.loc[method_order]
    rmse_pivot = rmse_pivot.loc[method_order]

    fig, axes = plt.subplots(1, 2, figsize=(18, max(8, 0.45 * len(method_order))))

    sns.heatmap(
        bias_pivot,
        annot=True,
        fmt=".3f",
        cmap="RdBu_r",
        center=0,
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "Mean Bias"},
        ax=axes[0],
    )
    axes[0].set_title("A) Bias Heatmap (True Slopes = 0.6, 1.0, 1.4)")
    axes[0].set_xlabel("True Slope (β)")
    axes[0].set_ylabel("Method")

    sns.heatmap(
        rmse_pivot,
        annot=True,
        fmt=".3f",
        cmap="YlOrRd",
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "RMSE"},
        ax=axes[1],
    )
    axes[1].set_title("B) RMSE Heatmap (True Slopes = 0.6, 1.0, 1.4)")
    axes[1].set_xlabel("True Slope (β)")
    axes[1].set_ylabel("Method")

    plt.tight_layout()
    out = OUT_DIR / "heatmap_bias_rmse_exact_slopes_single_habitat.png"
    plt.savefig(out, bbox_inches="tight")
    plt.show()
    print(f"Saved: {out}")


def plot_bias_rmse_heatmaps_exact_error_ratios(summary_df):
    plot_df = closest_exact_filter(summary_df, "error_ratio", ERROR_RATIOS_TARGET)

    bias_pivot = (
        plot_df.groupby(["method", "error_ratio"])["mean_bias"]
        .mean()
        .unstack()
        .reindex(columns=ERROR_RATIOS_TARGET)
    )

    rmse_pivot = (
        plot_df.groupby(["method", "error_ratio"])["rmse"]
        .mean()
        .unstack()
        .reindex(columns=ERROR_RATIOS_TARGET)
    )

    bias_pivot = bias_pivot.dropna(how="all", axis=0)
    rmse_pivot = rmse_pivot.dropna(how="all", axis=0)

    method_order = rmse_pivot.mean(axis=1).sort_values().index
    bias_pivot = bias_pivot.loc[method_order]
    rmse_pivot = rmse_pivot.loc[method_order]

    fig, axes = plt.subplots(1, 2, figsize=(18, max(8, 0.45 * len(method_order))))

    sns.heatmap(
        bias_pivot,
        annot=True,
        fmt=".3f",
        cmap="RdBu_r",
        center=0,
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "Mean Bias"},
        ax=axes[0],
    )
    axes[0].set_title("A) Bias Heatmap (Error Ratios = 0.5, 1.0, 2.0)")
    axes[0].set_xlabel("Error Ratio (σy / σx)")
    axes[0].set_ylabel("Method")

    sns.heatmap(
        rmse_pivot,
        annot=True,
        fmt=".3f",
        cmap="YlOrRd",
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "RMSE"},
        ax=axes[1],
    )
    axes[1].set_title("B) RMSE Heatmap (Error Ratios = 0.5, 1.0, 2.0)")
    axes[1].set_xlabel("Error Ratio (σy / σx)")
    axes[1].set_ylabel("Method")

    plt.tight_layout()
    out = OUT_DIR / "heatmap_bias_rmse_exact_error_ratios_single_habitat.png"
    plt.savefig(out, bbox_inches="tight")
    plt.show()
    print(f"Saved: {out}")


def plot_overall_bias_rmse_bars(summary_df):
    overall = (
        summary_df.groupby("method")
        .agg(mean_bias=("mean_bias", "mean"), rmse=("rmse", "mean"))
        .reset_index()
        .sort_values("rmse")
    )

    fig, axes = plt.subplots(1, 2, figsize=(18, max(8, 0.45 * len(overall))))

    bias_df = overall.sort_values("mean_bias")
    axes[0].barh(bias_df["method"], bias_df["mean_bias"], edgecolor="black", alpha=0.85)
    axes[0].axvline(0, color="black", linestyle="--", linewidth=1.5)
    axes[0].set_title("A) Overall Bias by Method")
    axes[0].set_xlabel("Mean Bias")

    rmse_df = overall.sort_values("rmse")
    axes[1].barh(rmse_df["method"], rmse_df["rmse"], edgecolor="black", alpha=0.85)
    axes[1].set_title("B) Overall RMSE by Method")
    axes[1].set_xlabel("RMSE")

    plt.tight_layout()
    out = OUT_DIR / "overall_bias_rmse_barplots_single_habitat.png"
    plt.savefig(out, bbox_inches="tight")
    plt.show()
    print(f"Saved: {out}")


def plot_rmse_vs_sample_size(summary_df, top_n=10):
    top_methods = (
        summary_df.groupby("method")["rmse"]
        .mean()
        .sort_values()
        .head(top_n)
        .index
        .tolist()
    )

    plot_df = summary_df[summary_df["method"].isin(top_methods)].copy()
    plot_df = closest_exact_filter(plot_df, "n", N_TARGET)

    fig, ax = plt.subplots(figsize=(10, 7))
    for method in top_methods:
        sub = plot_df[plot_df["method"] == method]
        line = sub.groupby("n")["rmse"].mean().reindex(N_TARGET)
        ax.plot(line.index, line.values, marker="o", linewidth=2, label=method)

    ax.set_title("RMSE vs Sample Size (Single Habitat)")
    ax.set_xlabel("Sample Size")
    ax.set_ylabel("RMSE")
    ax.set_xticks(N_TARGET)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=True)

    plt.tight_layout()
    out = OUT_DIR / "rmse_vs_sample_size_single_habitat.png"
    plt.savefig(out, bbox_inches="tight")
    plt.show()
    print(f"Saved: {out}")


def build_decision_tree_best_method(summary_df):
    scenario_cols = ["beta_true", "error_ratio", "n", "sigma_x_me", "sigma_y_me"]
    scenario_best = summary_df.loc[summary_df.groupby(scenario_cols)["rmse"].idxmin()].copy().reset_index(drop=True)

    X = scenario_best[["beta_true", "error_ratio", "n"]].copy()
    y = scenario_best["method"].copy()

    keep = X.notna().all(axis=1) & y.notna()
    X = X.loc[keep]
    y = y.loc[keep]

    if X.empty:
        print("No data available for decision tree.")
        return

    le = LabelEncoder()
    y_enc = le.fit_transform(y)

    clf = DecisionTreeClassifier(
        criterion="gini",
        max_depth=4,
        min_samples_leaf=2,
        random_state=42,
    )
    clf.fit(X, y_enc)

    fig, ax = plt.subplots(figsize=(28, 16))
    plot_tree(
        clf,
        feature_names=list(X.columns),
        class_names=list(le.classes_),
        filled=True,
        rounded=True,
        fontsize=14,
        ax=ax,
    )
    ax.set_title("Decision Tree for Choosing the Best Method (Single Habitat)")
    plt.tight_layout()

    out = OUT_DIR / "decision_tree_best_method_single_habitat.png"
    plt.savefig(out, bbox_inches="tight")
    plt.show()
    print(f"Saved: {out}")

    rules = export_text(clf, feature_names=list(X.columns))
    txt_out = OUT_DIR / "decision_tree_rules_single_habitat.txt"
    with open(txt_out, "w", encoding="utf-8") as f:
        f.write("Decision Tree Rules for Best Method (Single Habitat)\n")
        f.write("===================================================\n\n")
        f.write(rules)
        f.write("\n\nClass mapping:\n")
        for i, c in enumerate(le.classes_):
            f.write(f"{i}: {c}\n")
    print(f"Saved: {txt_out}")


def make_hierarchical_scenario_label(df):
    df = df.copy()
    df["beta_true_str"] = df["beta_true"].map(lambda x: f"β={x:g}")
    df["n_str"] = df["n"].map(lambda x: f"n={int(x)}")
    df["error_ratio_str"] = df["error_ratio"].map(lambda x: f"r={x:g}")
    return df


def plot_structured_hierarchical_heatmap(summary_df, value_col="rmse", suffix="rmse"):
    df = summary_df.copy()
    df = closest_exact_filter(df, "beta_true", TRUE_SLOPES_TARGET)
    df = closest_exact_filter(df, "n", N_TARGET)
    df = closest_exact_filter(df, "error_ratio", ERROR_RATIOS_TARGET)

    if df.empty:
        print(f"No data available for structured hierarchical heatmap ({value_col}).")
        return

    df = (
        df.groupby(["method", "beta_true", "n", "error_ratio"])[value_col]
        .mean()
        .reset_index()
    )

    df = df[~df["method"].isin(EXCLUDED_METHOD_NAMES)].copy()
    df = make_hierarchical_scenario_label(df)

    pivot = df.pivot_table(
        index="method",
        columns=["beta_true_str", "n_str", "error_ratio_str"],
        values=value_col,
        aggfunc="mean",
    )

    desired_cols = []
    for b in TRUE_SLOPES_TARGET:
        for n in N_TARGET:
            for er in ERROR_RATIOS_TARGET:
                desired_cols.append((f"β={b:g}", f"n={n}", f"r={er:g}"))

    existing_cols = [c for c in desired_cols if c in pivot.columns]
    pivot = pivot.reindex(columns=existing_cols)
    pivot = pivot.dropna(how="all", axis=0)

    method_order = pivot.mean(axis=1).sort_values().index
    pivot = pivot.loc[method_order]

    fig_w = max(18, 0.9 * len(pivot.columns))
    fig_h = max(8, 0.42 * len(pivot.index))

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    cmap = "YlOrRd" if value_col == "rmse" else "RdBu_r"
    center = 0 if value_col == "mean_bias" else None

    sns.heatmap(
        pivot,
        annot=True,
        fmt=".3f",
        cmap=cmap,
        center=center,
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"label": "RMSE" if value_col == "rmse" else "Mean Bias"},
        ax=ax,
    )

    title_metric = "RMSE" if value_col == "rmse" else "Bias"
    ax.set_title(
        f"Structured Hierarchical Heatmap ({title_metric})\n"
        f"β ∈ {TRUE_SLOPES_TARGET}, n ∈ {N_TARGET}, error ratio ∈ {ERROR_RATIOS_TARGET}"
    )
    ax.set_xlabel("Scenario Hierarchy: True Slope → Sample Size → Error Ratio")
    ax.set_ylabel("Method")

    plt.tight_layout()
    out = OUT_DIR / f"hierarchical_heatmap_{suffix}_single_habitat.png"
    plt.savefig(out, bbox_inches="tight")
    plt.show()
    print(f"Saved: {out}")


def create_summary_tables(summary_df):
    overall = (
        summary_df.groupby("method")
        .agg(
            mean_bias=("mean_bias", "mean"),
            rmse=("rmse", "mean"),
            mean_abs_bias=("mean_bias", lambda x: np.mean(np.abs(x))),
            n_scenarios=("rmse", "size"),
        )
        .sort_values("rmse")
        .reset_index()
    )

    by_slope = (
        summary_df.groupby(["method", "beta_true"])
        .agg(mean_bias=("mean_bias", "mean"), rmse=("rmse", "mean"))
        .reset_index()
    )

    by_error_ratio = (
        summary_df.groupby(["method", "error_ratio"])
        .agg(mean_bias=("mean_bias", "mean"), rmse=("rmse", "mean"))
        .reset_index()
    )

    by_n = (
        summary_df.groupby(["method", "n"])
        .agg(mean_bias=("mean_bias", "mean"), rmse=("rmse", "mean"))
        .reset_index()
    )

    return overall, by_slope, by_error_ratio, by_n

# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 90)
    print("SINGLE-HABITAT META-ANALYSIS PIPELINE")
    print("=" * 90)
    print(f"Output folder: {OUT_DIR}")

    df = load_all_single_habitat_data(BASE_DIR)
    print(f"Raw shape: {df.shape}")

    df = prepare_dataframe(df)
    print("Prepared dataframe.")
    print("Unique true slopes found:", sorted(df["beta_true"].dropna().unique().tolist()))
    print("Unique sample sizes found:", sorted(df["n"].dropna().unique().tolist()))
    print("Unique error ratios found:", sorted(np.round(df["error_ratio"].dropna().unique(), 8).tolist()))

    method_cols = find_method_columns(df)

    if not method_cols:
        raise ValueError("No beta_* method columns found in the data.")

    print("\nMethods found before final safeguard:")
    for m in method_cols:
        print(f"  - {m}")

    method_cols = [m for m in method_cols if m not in EXCLUDED_METHOD_COLUMNS]

    print(f"\nMethods used in analysis ({len(method_cols)}):")
    for m in method_cols:
        print(f"  - {m} -> {build_method_name(m)}")

    long_df = compute_long_metrics(df, method_cols)
    save_table(long_df, "long_metrics_rowwise_single_habitat")

    summary_df = summarize_metrics(long_df)

    # Final safeguard before tables/plots
    summary_df = summary_df[~summary_df["method_col"].isin(EXCLUDED_METHOD_COLUMNS)].copy()
    summary_df = summary_df[~summary_df["method"].isin(EXCLUDED_METHOD_NAMES)].copy()

    save_table(summary_df, "scenario_summary_metrics_single_habitat")

    overall, by_slope, by_error_ratio, by_n = create_summary_tables(summary_df)
    save_table(overall, "overall_method_summary_single_habitat")
    save_table(by_slope, "summary_by_slope_single_habitat")
    save_table(by_error_ratio, "summary_by_error_ratio_single_habitat")
    save_table(by_n, "summary_by_n_single_habitat")

    save_excel({
        "overall_summary": overall,
        "by_slope": by_slope,
        "by_error_ratio": by_error_ratio,
        "by_n": by_n,
        "scenario_summary": summary_df,
    })

    plot_bias_rmse_heatmaps_exact_slopes(summary_df)
    plot_bias_rmse_heatmaps_exact_error_ratios(summary_df)
    plot_overall_bias_rmse_bars(summary_df)
    plot_rmse_vs_sample_size(summary_df, top_n=10)
    build_decision_tree_best_method(summary_df)
    plot_structured_hierarchical_heatmap(summary_df, value_col="rmse", suffix="rmse")
    plot_structured_hierarchical_heatmap(summary_df, value_col="mean_bias", suffix="bias")

    print("\nDone.")
    print(f"All outputs saved in:\n{OUT_DIR}")


if __name__ == "__main__":
    main()
