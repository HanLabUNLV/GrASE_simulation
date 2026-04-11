#!/usr/bin/env python3
"""
Visualize GrASE evaluation results.

Three comparison groups:
  1. bipartition  -- all bipartition models and combination methods
  2. n_choose_2   -- all n_choose_2 models and combination methods
  3. cross        -- best bipartition vs best n_choose_2 vs multinomial
                     (best = highest micro_f1 at --padj for sim_type ALL)

Restricted evaluation plots are always generated alongside non-restricted,
in the same output subdirectory with a _restricted filename suffix.

Output layout:
  <out>/bipartition/  01_pr_curves.png  01_pr_curves_restricted.png  ...
  <out>/n_choose_2/   ...
  <out>/cross/        ...

Usage:
  python visualize_eval.py [--results-dir ~/GrASE_simulation/results] [--out plots/]
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns

# ---------------------------------------------------------------------------
# Method registry
# ---------------------------------------------------------------------------

BIPARTITION_METHODS = {
    "fixedEB":             {"dir": "eval_bipartition_fixedEB",             "prefix": "grase"},
    "fixedEB_mincomb":     {"dir": "eval_bipartition_fixedEB_mincomb",     "prefix": "grase"},
    "fixedEB_fishercomb":  {"dir": "eval_bipartition_fixedEB_fishercomb",  "prefix": "grase"},
    "prior":               {"dir": "eval_bipartition_prior",               "prefix": "grase"},
    "prior_mincomb":       {"dir": "eval_bipartition_prior_mincomb",       "prefix": "grase"},
    "prior_fishercomb":    {"dir": "eval_bipartition_prior_fishercomb",    "prefix": "grase"},
    "no_prior":            {"dir": "eval_bipartition_no_prior",            "prefix": "grase"},
    "no_prior_mincomb":    {"dir": "eval_bipartition_no_prior_mincomb",    "prefix": "grase"},
    "no_prior_fishercomb": {"dir": "eval_bipartition_no_prior_fishercomb", "prefix": "grase"},
    "wilcoxon":            {"dir": "eval_bipartition_wilcoxon",            "prefix": "grase"},
    "wilcoxon_mincomb":    {"dir": "eval_bipartition_wilcoxon_mincomb",    "prefix": "grase"},
    "wilcoxon_fishercomb": {"dir": "eval_bipartition_wilcoxon_fishercomb", "prefix": "grase"},
}

NC2_METHODS = {
    "fixedEB":             {"dir": "eval_n_choose_2_fixedEB",             "prefix": "grase"},
    "fixedEB_mincomb":     {"dir": "eval_n_choose_2_fixedEB_mincomb",     "prefix": "grase"},
    "fixedEB_fishercomb":  {"dir": "eval_n_choose_2_fixedEB_fishercomb",  "prefix": "grase"},
    "prior":               {"dir": "eval_n_choose_2_prior",               "prefix": "grase"},
    "prior_mincomb":       {"dir": "eval_n_choose_2_prior_mincomb",       "prefix": "grase"},
    "prior_fishercomb":    {"dir": "eval_n_choose_2_prior_fishercomb",    "prefix": "grase"},
    "no_prior":            {"dir": "eval_n_choose_2_no_prior",            "prefix": "grase"},
    "no_prior_mincomb":    {"dir": "eval_n_choose_2_no_prior_mincomb",    "prefix": "grase"},
    "no_prior_fishercomb": {"dir": "eval_n_choose_2_no_prior_fishercomb", "prefix": "grase"},
    "wilcoxon":            {"dir": "eval_n_choose_2_wilcoxon",            "prefix": "grase"},
    "wilcoxon_mincomb":    {"dir": "eval_n_choose_2_wilcoxon_mincomb",    "prefix": "grase"},
    "wilcoxon_fishercomb": {"dir": "eval_n_choose_2_wilcoxon_fishercomb", "prefix": "grase"},
}

MULTINOMIAL_METHODS = {
    "plugin_dm_EB": {"dir": "eval_multinomial_plugin_dm", "prefix": "grase"},
}

TOOLS_METHODS = {
    "rmats":             {"dir": "eval_saturn_rmats", "prefix": "rmats"},
    "dexseq":            {"dir": "eval_dexseq_rmats", "prefix": "dexseq"},
    "saturn":            {"dir": "eval_saturn_rmats", "prefix": "saturn"},
    "saturn_regularFDR": {"dir": "eval_saturn_rmats", "prefix": "saturn_regularFDR"},
}

# All GrASE methods with namespaced keys for best-grase selection
ALL_GRASE_METHODS = {}
for _k, _v in BIPARTITION_METHODS.items():
    ALL_GRASE_METHODS[f"bipartition.{_k}"] = _v
for _k, _v in NC2_METHODS.items():
    ALL_GRASE_METHODS[f"n_choose_2.{_k}"] = _v
for _k, _v in MULTINOMIAL_METHODS.items():
    ALL_GRASE_METHODS[f"multinomial.{_k}"] = _v

PADJ_THRESHOLDS = [0.01, 0.05, 0.10, 0.20]
SIM_TYPES_ORDERED = ["ALL", "Background", "DGE", "DTE", "DTU"]

MODEL_LABELS = {
    "fixedEB":             "EBfixed",
    "fixedEB_mincomb":     "EBfixed (mincomb)",
    "fixedEB_fishercomb":  "EBfixed (fishercomb)",
    "prior":               "EBprior",
    "prior_mincomb":       "EBprior (mincomb)",
    "prior_fishercomb":    "EBprior (fishercomb)",
    "no_prior":            "MLE",
    "no_prior_mincomb":    "MLE (mincomb)",
    "no_prior_fishercomb": "MLE (fishercomb)",
    "wilcoxon":            "wilcoxon",
    "wilcoxon_mincomb":    "wilcoxon (mincomb)",
    "wilcoxon_fishercomb": "wilcoxon (fishercomb)",
    "plugin_dm_EB":        "Multinomial (plugin_dm_EB)",
    "rmats":               "rMATS",
    "dexseq":              "DEXSeq",
    "saturn":              "Saturn (empirical FDR)",
    "saturn_regularFDR":   "Saturn (regular FDR)",
}

# Labels for namespaced GrASE keys used in tools comparison
def _grase_label(ns_key):
    """bipartition.fixedEB -> 'GrASE bipartition fixedEB'"""
    parts = ns_key.split(".", 1)
    group = parts[0]
    mkey  = parts[1] if len(parts) > 1 else ""
    base  = MODEL_LABELS.get(mkey, mkey)
    return f"GrASE {group} {base}"

# Within-group palette: each base model gets a color family (light/medium/dark)
# fixedEB = blue, prior = red/orange, no_prior = green, wilcoxon = purple
WITHIN_PALETTE = {
    "fixedEB":             "#1f77b4",  # dark blue
    "fixedEB_mincomb":     "#74b9d9",  # mid blue
    "fixedEB_fishercomb":  "#aec6e8",  # light blue
    "prior":               "#d62728",  # dark red
    "prior_mincomb":       "#eb8a5a",  # mid orange-red
    "prior_fishercomb":    "#f5a97a",  # light orange
    "no_prior":            "#2ca02c",  # dark green
    "no_prior_mincomb":    "#74c274",  # mid green
    "no_prior_fishercomb": "#a8d5a2",  # light green
    "wilcoxon":            "#9467bd",  # dark purple
    "wilcoxon_mincomb":    "#ad8dcc",  # mid purple
    "wilcoxon_fishercomb": "#c9b3d9",  # light purple
}

# Fixed colors for the cross-group comparison
CROSS_PALETTE = {
    "bipartition": "#1f77b4",
    "n_choose_2":  "#d62728",
    "multinomial": "#2ca02c",
}

TOOLS_PALETTE = {
    "rmats":             "#e377c2",  # pink
    "dexseq":            "#9467bd",  # purple
    "saturn":            "#8c564b",  # brown
    "saturn_regularFDR": "#c49c94",  # light brown
}

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _summary_path(results_dir, methods_dict, method, restricted):
    cfg  = methods_dict[method]
    tag  = "restricted_summary" if restricted else "summary"
    return results_dir / cfg["dir"] / f"{cfg['prefix']}_{tag}_by_simtype.txt"


def _per_gene_path(results_dir, methods_dict, method, padj, restricted):
    cfg  = methods_dict[method]
    tag  = f"restricted_per_gene_padj{padj:.2f}" if restricted else f"per_gene_padj{padj:.2f}"
    return results_dir / cfg["dir"] / f"{cfg['prefix']}_{tag}.txt"


def load_summaries(results_dir, methods_dict, restricted=False):
    frames = []
    for method in methods_dict:
        p = _summary_path(results_dir, methods_dict, method, restricted)
        if not p.exists():
            print(f"  [warn] missing: {p}")
            continue
        df = pd.read_csv(p, sep="\t")
        df["method"] = method
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    for col in ["micro_precision", "micro_recall", "micro_f1",
                "macro_precision", "macro_recall", "macro_f1"]:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def load_per_gene(results_dir, methods_dict, restricted=False):
    frames = []
    for method in methods_dict:
        for padj in PADJ_THRESHOLDS:
            p = _per_gene_path(results_dir, methods_dict, method, padj, restricted)
            if not p.exists():
                continue
            df = pd.read_csv(p, sep="\t")
            df["method"] = method
            frames.append(df)
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    for col in ["precision", "recall", "f1"]:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def best_method(summary, padj=0.05):
    """Return method with highest micro_f1 at padj for sim_type ALL."""
    sub = summary[(summary["sim_type"] == "ALL") & (summary["padj_thr"] == padj)]
    if sub.empty:
        return None
    return sub.loc[sub["micro_f1"].idxmax(), "method"]

# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def _save(fig, out_dir, name):
    path = out_dir / f"{name}.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved: {path}")


def _legend_handles(methods, palette, labels):
    return [
        plt.Line2D([0], [0], color=palette[m], lw=2.5, label=labels.get(m, m))
        for m in methods if m in palette
    ]

# ---------------------------------------------------------------------------
# Figure functions
# ---------------------------------------------------------------------------

def fig_pr_curves(summary, methods, palette, labels, out_dir, title_prefix, suffix=""):
    sim_types = [s for s in SIM_TYPES_ORDERED
                 if s in summary["sim_type"].unique() and s not in ("DGE", "Background")]
    ncols = len(sim_types)
    fig, axes = plt.subplots(1, ncols, figsize=(4 * ncols, 4))
    if ncols == 1:
        axes = [axes]

    for ax, stype in zip(axes, sim_types):
        # find the method with the lowest mean recall to annotate thresholds on
        mean_prec = {}
        for method in methods:
            sub = summary[(summary["sim_type"] == stype) & (summary["method"] == method)]
            if not sub.empty:
                mean_prec[method] = sub["micro_precision"].mean()
        label_method = min(mean_prec, key=mean_prec.get) if mean_prec else None

        for method in methods:
            sub = (summary[(summary["sim_type"] == stype) &
                           (summary["method"] == method)]
                   .sort_values("padj_thr"))
            if sub.empty:
                continue
            prec = sub["micro_precision"].values
            rec  = sub["micro_recall"].values
            ax.plot(rec, prec, "o-", color=palette[method], alpha=0.6,
                    label=labels.get(method, method), lw=2, ms=6)
            if method == label_method:
                for thr, r, p in zip(sub["padj_thr"], rec, prec):
                    ax.annotate(f"{thr}", (r, p), textcoords="offset points",
                                xytext=(4, -8), fontsize=6, color="black")
        ax.set_title(stype, fontsize=12, fontweight="bold")
        ax.set_xlabel("Recall", fontsize=10)
        ax.set_ylabel("Precision", fontsize=10)
        # zoom to data range with a small margin
        all_rec  = summary[summary["sim_type"] == stype]["micro_recall"].dropna()
        all_prec = summary[summary["sim_type"] == stype]["micro_precision"].dropna()
        if len(all_rec) > 0 and len(all_prec) > 0:
            pad = 0.03
            ax.set_xlim(max(0, all_rec.min()  - pad), min(1, all_rec.max()  + pad))
            ax.set_ylim(max(0, all_prec.min() - pad), min(1, all_prec.max() + pad))
        else:
            ax.set_xlim(-0.02, 1.02)
            ax.set_ylim(-0.02, 1.02)
        ax.plot([0, 1], [0, 1], "k--", lw=0.8, alpha=0.4)
        ax.grid(True, alpha=0.3)

    ncol_leg = min(4, len(methods))
    fig.legend(handles=_legend_handles(methods, palette, labels),
               loc="lower center", ncol=ncol_leg,
               fontsize=8, bbox_to_anchor=(0.5, -0.12))
    rest_tag = " (restricted)" if suffix else ""
    fig.suptitle(f"{title_prefix} -- Precision-Recall{rest_tag}",
                 fontsize=13, fontweight="bold", y=1.02)
    fig.tight_layout()
    _save(fig, out_dir, f"01_pr_curves{suffix}")


def _group_by_base(methods):
    """Group methods by base model (strip _mincomb / _fishercomb suffixes)."""
    groups = {}
    for m in methods:
        base = m.replace("_fishercomb", "").replace("_mincomb", "")
        groups.setdefault(base, []).append(m)
    return groups


def fig_pr_curves_panels(summary, methods, palette, labels, out_dir, title_prefix, suffix=""):
    sim_types = [s for s in SIM_TYPES_ORDERED
                 if s in summary["sim_type"].unique() and s not in ("DGE", "Background")]

    # group by combination variant: plain, mincomb, fishercomb
    comb_groups = {
        "plain":      [m for m in methods if not m.endswith("_mincomb")
                                         and not m.endswith("_fishercomb")],
        "mincomb":    [m for m in methods if m.endswith("_mincomb")],
        "fishercomb": [m for m in methods if m.endswith("_fishercomb")],
    }
    comb_names = [c for c in comb_groups if comb_groups[c]]

    # all variants of a base model share the plain color
    panel_palette = {}
    for m in methods:
        base = m.replace("_fishercomb", "").replace("_mincomb", "")
        panel_palette[m] = palette[base]

    nrows = len(comb_names)
    ncols = len(sim_types)
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(4 * ncols, 3.5 * nrows),
                             squeeze=False)

    # compute per-sim_type axis limits so each column zooms to its own data range
    pad = 0.03
    stype_lims = {}
    for stype in sim_types:
        sub = summary[summary["sim_type"] == stype]
        r = sub["micro_recall"].dropna()
        p = sub["micro_precision"].dropna()
        if len(r) > 0 and len(p) > 0:
            stype_lims[stype] = (
                (max(0, r.min() - pad), min(1, r.max() + pad)),
                (max(0, p.min() - pad), min(1, p.max() + pad)),
            )
        else:
            stype_lims[stype] = ((-0.02, 1.02), (-0.02, 1.02))

    for ri, comb in enumerate(comb_names):
        comb_methods = comb_groups[comb]
        for ci, stype in enumerate(sim_types):
            ax = axes[ri][ci]

            mean_prec = {}
            for method in comb_methods:
                sub = summary[(summary["sim_type"] == stype) &
                              (summary["method"] == method)]
                if not sub.empty:
                    mean_prec[method] = sub["micro_precision"].mean()
            label_method = min(mean_prec, key=mean_prec.get) if mean_prec else None

            for method in comb_methods:
                sub = (summary[(summary["sim_type"] == stype) &
                               (summary["method"] == method)]
                       .sort_values("padj_thr"))
                if sub.empty:
                    continue
                prec = sub["micro_precision"].values
                rec  = sub["micro_recall"].values
                ax.plot(rec, prec, "o-", color=panel_palette[method], alpha=0.7,
                        label=labels.get(method, method), lw=2, ms=6)
                if method == label_method:
                    for thr, r, p in zip(sub["padj_thr"], rec, prec):
                        ax.annotate(f"{thr}", (r, p), textcoords="offset points",
                                    xytext=(4, -8), fontsize=6, color="black")

            xlim, ylim = stype_lims[stype]
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            ax.plot([0, 1], [0, 1], "k--", lw=0.8, alpha=0.4)
            ax.grid(True, alpha=0.3)

            if ri == 0:
                ax.set_title(stype, fontsize=11, fontweight="bold")
            if ci == 0:
                ax.set_ylabel(f"{comb}\nPrecision", fontsize=10, fontweight="bold")
            ax.set_xlabel("Recall" if ri == nrows - 1 else "", fontsize=9)

            ax.legend(fontsize=7, loc="best",
                      handles=_legend_handles(comb_methods, panel_palette, labels))

    rest_tag = " (restricted)" if suffix else ""
    fig.suptitle(f"{title_prefix} -- Precision-Recall by combination method{rest_tag}",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    _save(fig, out_dir, f"01b_pr_curves_panels{suffix}")


def fig_f1_heatmap(summary, methods, labels, out_dir, title_prefix, suffix=""):
    sim_types = [s for s in SIM_TYPES_ORDERED if s in summary["sim_type"].unique()]
    nrows = len(PADJ_THRESHOLDS)
    fig, axes = plt.subplots(nrows, 1, figsize=(len(methods) * 1.5 + 2, nrows * 2.5))
    if nrows == 1:
        axes = [axes]

    for ax, thr in zip(axes, PADJ_THRESHOLDS):
        sub = summary[summary["padj_thr"] == thr]
        pivot = (sub.pivot_table(index="method", columns="sim_type",
                                 values="micro_f1", aggfunc="first")
                    .reindex(index=methods, columns=sim_types))
        pivot.index = [labels.get(m, m) for m in pivot.index]
        sns.heatmap(pivot, ax=ax, annot=True, fmt=".3f", cmap="YlOrRd",
                    vmin=0, vmax=1, linewidths=0.5, cbar=True,
                    cbar_kws={"shrink": 0.8})
        ax.set_title(f"padj <= {thr}", fontsize=10)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(axis="x", labelsize=9)
        ax.tick_params(axis="y", labelsize=9, rotation=0)

    rest_tag = " (restricted)" if suffix else ""
    fig.suptitle(f"{title_prefix} -- Micro F1{rest_tag}",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    _save(fig, out_dir, f"02_f1_heatmap{suffix}")


def fig_metric_bars(summary, methods, palette, labels, out_dir,
                    title_prefix, padj=0.05, suffix=""):
    sim_types = [s for s in SIM_TYPES_ORDERED
                 if s in summary["sim_type"].unique() and s != "DGE"]
    metrics      = ["micro_precision", "micro_recall", "micro_f1"]
    metric_labels = ["Precision", "Recall", "F1"]

    nrows = len(sim_types)
    ncols = len(metrics)
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * max(5, len(methods) * 0.6),
                                                     nrows * 3))
    sub = summary[summary["padj_thr"] == padj]
    x   = np.arange(len(methods))
    width = 0.6

    for ri, stype in enumerate(sim_types):
        row_data = sub[sub["sim_type"] == stype]
        for ci, (metric, mlabel) in enumerate(zip(metrics, metric_labels)):
            ax = axes[ri][ci] if nrows > 1 else axes[ci]
            vals = [row_data.loc[row_data["method"] == m, metric].values for m in methods]
            vals = [v[0] if len(v) > 0 and not np.isnan(v[0]) else 0.0 for v in vals]
            bars = ax.bar(x, vals, width, color=[palette[m] for m in methods],
                          edgecolor="white", lw=0.5)
            ax.set_xticks(x)
            ax.set_xticklabels([labels.get(m, m) for m in methods],
                               rotation=35, ha="right", fontsize=8)
            ax.set_ylim(0, 1.05)
            ax.set_ylabel(mlabel if ci == 0 else "", fontsize=10)
            ax.set_title(f"{stype} -- {mlabel}", fontsize=10)
            ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
            ax.grid(axis="y", alpha=0.3)
            for bar, val in zip(bars, vals):
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                        f"{val:.3f}", ha="center", va="bottom", fontsize=7)

    rest_tag = " (restricted)" if suffix else ""
    fig.suptitle(f"{title_prefix} -- Metrics at padj <= {padj}{rest_tag}",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    _save(fig, out_dir, f"03_metric_bars_padj{padj}{suffix}")


def fig_threshold_curves(summary, methods, palette, labels, out_dir,
                         title_prefix, suffix=""):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    metrics  = ["micro_precision", "micro_recall", "micro_f1"]
    mlabels  = ["Micro Precision", "Micro Recall", "Micro F1"]
    sub = summary[summary["sim_type"] == "ALL"].sort_values("padj_thr")

    for ax, metric, mlabel in zip(axes, metrics, mlabels):
        for method in methods:
            msub = sub[sub["method"] == method]
            ax.plot(msub["padj_thr"], msub[metric], "o-",
                    color=palette[method], label=labels.get(method, method),
                    lw=2, ms=6)
        ax.set_xlabel("padj threshold", fontsize=10)
        ax.set_ylabel(mlabel, fontsize=10)
        ax.set_title(mlabel, fontsize=11)
        ax.set_xlim(0, 0.22)
        ax.set_ylim(0, 1.0)
        ax.set_xticks(PADJ_THRESHOLDS)
        ax.grid(True, alpha=0.3)

    fig.legend(handles=_legend_handles(methods, palette, labels),
               loc="lower center", ncol=min(4, len(methods)),
               fontsize=8, bbox_to_anchor=(0.5, -0.12))
    rest_tag = " (restricted)" if suffix else ""
    fig.suptitle(f"{title_prefix} -- Metrics vs padj (ALL){rest_tag}",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    _save(fig, out_dir, f"04_threshold_curves{suffix}")


def fig_per_gene_boxplot(per_gene, methods, palette, labels, out_dir,
                         title_prefix, padj=0.05, suffix=""):
    sim_types = [s for s in SIM_TYPES_ORDERED
                 if s in per_gene["sim_type"].unique()
                 and s not in ("ALL", "Background", "DGE")]
    sub = per_gene[
        (per_gene["padj_thr"] == padj) &
        (per_gene["sim_type"].isin(sim_types)) &
        (per_gene["f1"].notna())
    ].copy()
    sub["method_label"] = sub["method"].map(labels)
    if sub.empty:
        print("  [skip] no per-gene data for boxplot")
        return

    order          = [labels.get(m, m) for m in methods]
    palette_mapped = {labels.get(m, m): palette[m] for m in methods}

    fig, axes = plt.subplots(1, len(sim_types),
                             figsize=(max(5, len(methods) * 0.6) * len(sim_types), 5))
    if len(sim_types) == 1:
        axes = [axes]

    for ax, stype in zip(axes, sim_types):
        sdata = sub[sub["sim_type"] == stype]
        sns.boxplot(data=sdata, x="method_label", y="f1", order=order,
                    hue="method_label", palette=palette_mapped, legend=False,
                    ax=ax, width=0.6,
                    flierprops=dict(marker=".", markersize=2, alpha=0.3))
        ax.set_title(stype, fontsize=12, fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel("Gene-level F1" if sim_types.index(stype) == 0 else "")
        ax.set_ylim(-0.05, 1.05)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order, rotation=35, ha="right", fontsize=8)
        ax.grid(axis="y", alpha=0.3)

    rest_tag = " (restricted)" if suffix else ""
    fig.suptitle(f"{title_prefix} -- Per-gene F1 at padj <= {padj}{rest_tag}",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    _save(fig, out_dir, f"05_per_gene_boxplot_padj{padj}{suffix}")


def fig_confusion_counts(summary, methods, palette, labels, out_dir,
                         title_prefix, padj=0.05, suffix=""):
    sim_types = [s for s in SIM_TYPES_ORDERED
                 if s in summary["sim_type"].unique()]
    sub   = summary[summary["padj_thr"] == padj]
    nrows = len(sim_types)
    fig, axes = plt.subplots(nrows, 1,
                             figsize=(max(8, len(methods) * 0.9), 3 * nrows))
    if nrows == 1:
        axes = [axes]

    x     = np.arange(len(methods))
    width = 0.6

    for ax, stype in zip(axes, sim_types):
        row = sub[sub["sim_type"] == stype]
        tp = [row.loc[row["method"] == m, "total_TP"].values for m in methods]
        fp = [row.loc[row["method"] == m, "total_FP"].values for m in methods]
        fn = [row.loc[row["method"] == m, "total_FN"].values for m in methods]
        tp = [v[0] if len(v) > 0 else 0 for v in tp]
        fp = [v[0] if len(v) > 0 else 0 for v in fp]
        fn = [v[0] if len(v) > 0 else 0 for v in fn]

        ax.bar(x, tp, width, label="TP", color="#2ca02c", edgecolor="white")
        ax.bar(x, fp, width, bottom=tp, label="FP", color="#d62728", edgecolor="white")
        ax.bar(x, fn, width, bottom=np.array(tp) + np.array(fp),
               label="FN", color="#ff7f0e", edgecolor="white")
        ax.set_xticks(x)
        ax.set_xticklabels([labels.get(m, m) for m in methods],
                           rotation=35, ha="right", fontsize=8)
        ax.set_ylabel("Count")
        ax.set_title(stype, fontsize=11)
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(axis="y", alpha=0.3)

    rest_tag = " (restricted)" if suffix else ""
    fig.suptitle(f"{title_prefix} -- TP/FP/FN at padj <= {padj}{rest_tag}",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    _save(fig, out_dir, f"06_confusion_counts_padj{padj}{suffix}")

# ---------------------------------------------------------------------------
# Run one comparison group (both restricted and non-restricted)
# ---------------------------------------------------------------------------

def run_group(results_dir, methods_dict, palette, labels, group_name,
              out_base, padj=0.05, panels=True):
    out_dir  = out_base / group_name
    out_dir.mkdir(parents=True, exist_ok=True)
    methods  = list(methods_dict.keys())
    title    = group_name.replace("_", " ")
    summary_full = pd.DataFrame()

    for restricted in (False, True):
        suffix = "_restricted" if restricted else ""
        tag    = "(restricted)" if restricted else "(full)"
        print(f"\n  Loading {group_name} {tag}...")
        summary  = load_summaries(results_dir, methods_dict, restricted=restricted)
        per_gene = load_per_gene(results_dir, methods_dict, restricted=restricted)
        if summary.empty:
            print(f"  [skip] no summary data for {group_name} {tag}")
            continue

        if not restricted:
            summary_full = summary

        print(f"  Generating {group_name} {tag} plots...")
        fig_pr_curves(summary, methods, palette, labels, out_dir, title, suffix=suffix)
        if panels:
            fig_pr_curves_panels(summary, methods, palette, labels, out_dir, title, suffix=suffix)
        fig_f1_heatmap(summary, methods, labels, out_dir, title, suffix=suffix)
        fig_metric_bars(summary, methods, palette, labels, out_dir, title,
                        padj=padj, suffix=suffix)
        fig_threshold_curves(summary, methods, palette, labels, out_dir, title,
                             suffix=suffix)
        if not per_gene.empty:
            fig_per_gene_boxplot(per_gene, methods, palette, labels, out_dir, title,
                                 padj=padj, suffix=suffix)
        fig_confusion_counts(summary, methods, palette, labels, out_dir, title,
                             padj=padj, suffix=suffix)
        # separate plots for raw vs combined-pval methods
        raw_methods  = [m for m in methods
                        if not any(c in m for c in ("mincomb", "fishercomb"))]
        comb_methods = [m for m in methods
                        if any(c in m for c in ("mincomb", "fishercomb"))]
        if raw_methods and comb_methods:
            fig_confusion_counts(summary, raw_methods, palette, labels, out_dir,
                                 title + " (raw)", padj=padj,
                                 suffix=suffix + "_raw")
            fig_confusion_counts(summary, comb_methods, palette, labels, out_dir,
                                 title + " (combined pval)", padj=padj,
                                 suffix=suffix + "_comb")

    return summary_full  # non-restricted, used by cross comparison

# ---------------------------------------------------------------------------
# Cross comparison: best bipartition vs best n_choose_2 vs multinomial
# ---------------------------------------------------------------------------

def run_cross(results_dir, summary_bip, summary_nc2, out_base, padj=0.05,
              user_best_bip=None, user_best_nc2=None):
    out_dir = out_base / "cross"
    out_dir.mkdir(parents=True, exist_ok=True)

    best_bip = user_best_bip if user_best_bip else (
        best_method(summary_bip, padj) if not summary_bip.empty else None)
    best_nc2 = user_best_nc2 if user_best_nc2 else (
        best_method(summary_nc2, padj) if not summary_nc2.empty else None)

    if best_bip is None or best_nc2 is None:
        print("  [skip] cannot determine best methods -- missing summary data "
              "(run bipartition and n_choose_2 groups first, or use "
              "--best-bipartition / --best-n-choose-2)")
        return

    src_bip = "user-specified" if user_best_bip else "auto-selected"
    src_nc2 = "user-specified" if user_best_nc2 else "auto-selected"

    print(f"\n  Best bipartition : {best_bip} ({src_bip})")
    print(f"  Best n_choose_2  : {best_nc2} ({src_nc2})")

    key_bip   = f"bipartition ({best_bip})"
    key_nc2   = f"n_choose_2 ({best_nc2})"
    key_multi = "multinomial (plugin_dm_EB)"

    cross_methods = {
        key_bip:   BIPARTITION_METHODS[best_bip],
        key_nc2:   NC2_METHODS[best_nc2],
        key_multi: MULTINOMIAL_METHODS["plugin_dm_EB"],
    }
    cross_palette = {
        key_bip:   CROSS_PALETTE["bipartition"],
        key_nc2:   CROSS_PALETTE["n_choose_2"],
        key_multi: CROSS_PALETTE["multinomial"],
    }
    cross_labels = {k: k for k in cross_methods}
    methods      = list(cross_methods.keys())

    for restricted in (False, True):
        suffix = "_restricted" if restricted else ""
        tag    = "(restricted)" if restricted else "(full)"
        print(f"\n  Loading cross {tag}...")
        summary  = load_summaries(results_dir, cross_methods, restricted=restricted)
        per_gene = load_per_gene(results_dir, cross_methods, restricted=restricted)
        if summary.empty:
            print(f"  [skip] no data for cross {tag}")
            continue

        print(f"  Generating cross {tag} plots...")
        fig_pr_curves(summary, methods, cross_palette, cross_labels,
                      out_dir, "Cross comparison", suffix=suffix)
        fig_f1_heatmap(summary, methods, cross_labels,
                       out_dir, "Cross comparison", suffix=suffix)
        fig_metric_bars(summary, methods, cross_palette, cross_labels,
                        out_dir, "Cross comparison", padj=padj, suffix=suffix)
        fig_threshold_curves(summary, methods, cross_palette, cross_labels,
                             out_dir, "Cross comparison", suffix=suffix)
        if not per_gene.empty:
            fig_per_gene_boxplot(per_gene, methods, cross_palette, cross_labels,
                                 out_dir, "Cross comparison", padj=padj, suffix=suffix)
        fig_confusion_counts(summary, methods, cross_palette, cross_labels,
                             out_dir, "Cross comparison", padj=padj, suffix=suffix)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description="Visualize GrASE evaluation results")
    parser.add_argument("--results-dir",
                        default=str(Path.home() / "GrASE_simulation/results"),
                        help="Results directory (default: ~/GrASE_simulation/results)")
    parser.add_argument("--out", default="plots",
                        help="Output directory for plots (default: ./plots)")
    parser.add_argument("--padj", type=float, default=0.05,
                        help="Primary padj threshold for bar/box plots and best-method "
                             "selection (default: 0.05)")
    parser.add_argument("--groups", nargs="+",
                        choices=["bipartition", "n_choose_2", "cross", "tools"],
                        default=["bipartition", "n_choose_2", "cross", "tools"],
                        help="Which comparison groups to run (default: all)")
    parser.add_argument("--best-bipartition", default=None,
                        choices=list(BIPARTITION_METHODS.keys()),
                        metavar="METHOD",
                        help="Bipartition method to use in cross comparison. "
                             "If omitted, auto-selected by highest micro_f1 at --padj. "
                             "Choices: " + ", ".join(BIPARTITION_METHODS.keys()))
    parser.add_argument("--best-n-choose-2", default=None,
                        choices=list(NC2_METHODS.keys()),
                        metavar="METHOD",
                        help="N-choose-2 method to use in cross comparison. "
                             "If omitted, auto-selected by highest micro_f1 at --padj. "
                             "Choices: " + ", ".join(NC2_METHODS.keys()))
    parser.add_argument("--best-grase", default=None,
                        choices=list(ALL_GRASE_METHODS.keys()),
                        metavar="GROUP.METHOD",
                        help="GrASE method to include in tools comparison, as "
                             "group.method (e.g. bipartition.fixedEB). "
                             "If omitted, auto-selected by highest micro_f1 at --padj "
                             "across all GrASE methods.")
    return parser.parse_args()


def main():
    args        = parse_args()
    results_dir = Path(args.results_dir).expanduser()
    out_base    = Path(args.out)
    out_base.mkdir(parents=True, exist_ok=True)

    print(f"Results dir : {results_dir}")
    print(f"Output dir  : {out_base}")
    print(f"Primary padj: {args.padj}")

    bip_palette   = {m: WITHIN_PALETTE[m] for m in BIPARTITION_METHODS}
    nc2_palette   = {m: WITHIN_PALETTE[m] for m in NC2_METHODS}
    bip_labels    = {m: MODEL_LABELS[m]   for m in BIPARTITION_METHODS}
    nc2_labels    = {m: MODEL_LABELS[m]   for m in NC2_METHODS}
    tools_labels  = {m: MODEL_LABELS[m]   for m in TOOLS_METHODS}

    groups = args.groups

    # 1. Bipartition comparison
    summary_bip = pd.DataFrame()
    if "bipartition" in groups:
        print("\n=== Bipartition ===")
        summary_bip = run_group(results_dir, BIPARTITION_METHODS, bip_palette, bip_labels,
                                "bipartition", out_base, padj=args.padj)

    # 2. N-choose-2 comparison
    summary_nc2 = pd.DataFrame()
    if "n_choose_2" in groups:
        print("\n=== N-choose-2 ===")
        summary_nc2 = run_group(results_dir, NC2_METHODS, nc2_palette, nc2_labels,
                                "n_choose_2", out_base, padj=args.padj)

    # 3. Tools comparison (rMATS / Saturn + best GrASE)
    if "tools" in groups:
        print("\n=== Tools (rMATS / Saturn + best GrASE) ===")

        # resolve best GrASE method
        if args.best_grase:
            best_grase_key = args.best_grase
            src = "user-specified"
        else:
            all_summary = load_summaries(results_dir, ALL_GRASE_METHODS, restricted=False)
            if all_summary.empty:
                best_grase_key = None
                print("  [warn] could not load any GrASE summaries for auto-selection")
            else:
                best_grase_key = best_method(all_summary, args.padj)
            src = "auto-selected"
        print(f"  Best GrASE : {best_grase_key} ({src})")

        # build augmented methods/palette/labels with best GrASE prepended
        tools_methods_aug  = dict(TOOLS_METHODS)
        tools_palette_aug  = dict(TOOLS_PALETTE)
        tools_labels_aug   = dict(tools_labels)
        if best_grase_key:
            tools_methods_aug  = {best_grase_key: ALL_GRASE_METHODS[best_grase_key],
                                   **TOOLS_METHODS}
            tools_palette_aug  = {best_grase_key: "#ff7f0e", **TOOLS_PALETTE}
            tools_labels_aug   = {best_grase_key: _grase_label(best_grase_key),
                                   **tools_labels}

        run_group(results_dir, tools_methods_aug, tools_palette_aug, tools_labels_aug,
                  "tools", out_base, padj=args.padj, panels=False)

    # 4. Cross comparison: best bipartition vs best n_choose_2 vs multinomial
    if "cross" in groups:
        print("\n=== Cross comparison ===")
        run_cross(results_dir, summary_bip, summary_nc2,
                  out_base=out_base, padj=args.padj,
                  user_best_bip=args.best_bipartition,
                  user_best_nc2=args.best_n_choose_2)

    print("\nDone.")


if __name__ == "__main__":
    main()
