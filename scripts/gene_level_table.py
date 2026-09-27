"""Gene-level metrics table, the gene-currency analogue of
results/eval_metricfair/transcript_level_metrics.txt.

A call implicates exactly ONE gene, so two columns of the transcript table
have no gene-level analogue and are deliberately omitted:
  - sharpness (implicated-set size) is always 1
  - prec_strict and prec_tolerant collapse into a single precision
    (a call either lands in a signal gene or a null gene)

A gene is CALLED if it carries >= 1 significant call; POSITIVE if DTE or DTU;
NEGATIVE if Background or DGE. Both universes are emitted:
  restricted -- denominator = signal genes the tool actually tested. Answers
                "when a tool can test a gene, how well does it do?" Recall is
                NOT comparable across tools (each has its own denominator).
  full       -- denominator = all signal genes (3,002; 1,501 per class).
                Structural non-coverage is charged as FN, so recall IS
                comparable across tools. This is the universe that credits
                reach, and is the one to lead with in cross-tool claims.
TP/FP/precision are identical between universes by construction; only FN and
recall differ.

Source: results/eval_metricfair/pr_three_levels.GT_rule.txt (level == "gene"),
written by scripts/pr_curves_three_levels_gtrule.R. Operating points match the
rest of the benchmark: padj/FDR < 0.01 for GrASE/DEXSeq/rMATS, P(|dPSI|>=C)
>= 0.95 for MAJIQ.

Usage: python3 scripts/gene_level_table.py
Out:   results/eval_metricfair/gene_level_table.{restricted,full}.txt
"""
import csv, os, collections

BASE = "/mnt/data1/home/mirahan/GrASE_simulation"
OUT = os.path.join(BASE, "results/eval_metricfair")
SWEEP = os.path.join(OUT, "pr_three_levels.GT_rule.txt")
TXTAB = os.path.join(OUT, "transcript_level_metrics.txt")

THR = {"GrASE": "0.01", "GrASE_merged_all": "0.01",
       "GrASE_merged_dpi0.1": "0.01", "GrASE_merged_dpi0.2": "0.01",
       "GrASE_merged_internal": "0.01",
       "GrASE_dpi0.1": "0.01", "GrASE_dpi0.2": "0.01",
       "GrASE_internal": "0.01", "MAJIQ_C0.10": "0.95", "MAJIQ_C0.20": "0.95",
       "DEXSeq": "0.01", "rMATS": "0.01",
       "rMATS_dpsi0.1": "0.01", "rMATS_dpsi0.2": "0.01"}
ORDER = ["GrASE", "GrASE_merged_all",
         "GrASE_dpi0.1", "GrASE_merged_dpi0.1",
         "GrASE_dpi0.2", "GrASE_merged_dpi0.2",
         "GrASE_internal", "GrASE_merged_internal",
         "MAJIQ_C0.10", "MAJIQ_C0.20", "DEXSeq",
         "rMATS", "rMATS_dpsi0.1", "rMATS_dpsi0.2"]
CATS = ["ALL", "DTE", "DTU", "Background", "DGE"]

sweep = list(csv.DictReader(open(SWEEP), delimiter="\t"))
# n_sig_calls comes from the transcript table (same calls, already tabulated)
calls = {}
for r in csv.DictReader(open(TXTAB), delimiter="\t"):
    if r["category"] == "ALL":
        calls[r["tool"]] = r["n_sig_calls"]

idx = collections.defaultdict(dict)
for r in sweep:
    if r["level"] != "gene":
        continue
    if r["tool"] in THR and r["thr"] == THR[r["tool"]]:
        idx[(r["category"], r["universe"])][r["tool"]] = r

HDR = ["category", "tool", "thr", "n_sig_calls", "n_pos_genes", "n_neg_genes",
       "TP", "FP", "FN", "TN", "recall", "precision", "f1"]

for uni in ("restricted", "full"):
    path = os.path.join(OUT, f"gene_level_table.{uni}.txt")
    with open(path, "w") as fh:
        fh.write("\t".join(HDR) + "\n")
        for cat in CATS:
            for t in ORDER:
                r = idx.get((cat, uni), {}).get(t)
                if r is None:
                    continue
                P, R = float(r["precision"]), float(r["recall"])
                f1 = 0.0 if P + R == 0 else 2 * P * R / (P + R)
                fh.write("\t".join([
                    cat, t, THR[t], calls.get(t, "NA"),
                    r["n_pos"], r["n_neg"], r["TP"], r["FP"], r["FN"], r["TN"],
                    f"{R:.4f}", f"{P:.4f}", f"{f1:.4f}"]) + "\n")
    print(f"wrote {path}")
