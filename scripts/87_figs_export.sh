#!/bin/bash
# Export every figure PDF to the two formats the submission needs.
#
#   EPS  for the journal.
#   SVG  for PowerPoint. PowerPoint CANNOT import EPS -- Office disabled it in a
#        2017 security update -- but Office 2016 and later insert SVG as true
#        vector, so it scales in a slide with no quality loss.
#   PNG  600 dpi, for anywhere vector is awkward: the dispersion scatters carry
#        29k points each and their SVG runs to 13-15 MB, which PowerPoint opens
#        but redraws slowly.
#
# PDF stays the master and both are derived from it. EPS has no transparency and
# R's postscript() device cannot do alpha at all, so drawing EPS directly would
# silently drop the shading in the dispersion scatters; pdftops flattens it
# instead. The plotting scripts now write PDF only.
#
# Usage: bash scripts/87_figs_export.sh
set -eu
cd /mnt/data1/home/mirahan/GrASE_simulation
mkdir -p figs/eps figs/svg figs/png600

export_one () {
  local pdf="$1" base
  base=$(basename "${pdf%.pdf}")
  pdftops   -eps "$pdf" "figs/eps/$base.eps"
  pdftocairo -svg "$pdf" "figs/svg/$base.svg"
  # -singlefile keeps the name as given; pdftocairo would otherwise append -1.
  pdftocairo -png -r 600 -singlefile "$pdf" "figs/png600/$base"
  printf "  %-40s eps %-6s svg %-6s png %s\n" "$base" \
         "$(stat -c %s "figs/eps/$base.eps" | numfmt --to=iec)" \
         "$(stat -c %s "figs/svg/$base.svg" | numfmt --to=iec)" \
         "$(stat -c %s "figs/png600/$base.png" | numfmt --to=iec)"
}

echo "=== main panels ==="
for f in plots/pr_three_levels.GT_rule*.pdf \
         plots/confusion_three_levels.GT_rule*.pdf \
         plots/pr_models.GT_rule.*.pdf plots/delta_sweep.GT_rule.*.pdf \
         plots/roc_partial.GT_rule.*.pdf \
         plots/pr_dte_dtu.*.pdf plots/pr_all.*.pdf plots/roc_dte_dtu.*.pdf \
         plots/cross/01_pr_curves*.pdf plots/cross/06_confusion_counts_padj0.01.pdf \
         plots/dte_fp_example_*.pdf; do
  [ -f "$f" ] && export_one "$f"
done

echo "=== supplementary dispersion and effect-size panels ==="
for f in scripts/plots/model_comparison/internalp1.pdf \
         scripts/plots/model_comparison/internalp2.pdf \
         scripts/plots/model_comparison/internalp7.pdf; do
  [ -f "$f" ] && export_one "$f"
done

echo "=== AIF1 splicing graph and transcript structure (Figure 1) ==="
# Written by ~/GrASE/Rpkg/scripts/plot_aif1_fig1.R (grase style_and_plot / plottx).
for f in figs/ENSG00000204472.13.pdf figs/ENSG00000204472.13.tx.pdf; do
  [ -f "$f" ] && export_one "$f" \
    || echo "  MISSING $f -- run Rscript ~/GrASE/Rpkg/scripts/plot_aif1_fig1.R"
done

echo "=== DICE gene case-study panels ==="
# Written by ~/DICE/regen_fig9_fig10.sh (plot_test_results.py). They live in the
# DICE tree because that is where their inputs are; only the exports land here.
for f in ~/DICE/fig_merged/AIF1_ev11.pdf ~/DICE/fig_merged/ARHGAP9_ev73.pdf; do
  [ -f "$f" ] && export_one "$f" || echo "  MISSING $f -- run ~/DICE/regen_fig9_fig10.sh"
done

echo "=== hand-built SVG figures (graphical abstract, toy graph, subgraph) ==="
# These are written as SVG directly, not through a plotting device, so there is
# no PDF master. cairosvg renders them to true vector EPS and PDF from the same
# source; rewriting the generators to emit PostScript natively would buy nothing.
CPY=/mnt/data1/home/mirahan/miniconda3/envs/pdf/bin/python
for f in figs/graphical_abstract.svg figs/toy_grase_graph.svg figs/aif1_subgraph_R14.svg; do
  [ -f "$f" ] || continue
  b=$(basename "${f%.svg}")
  "$CPY" -c "import cairosvg
cairosvg.svg2eps(url='$f', write_to='figs/eps/$b.eps')
cairosvg.svg2pdf(url='$f', write_to='figs/$b.pdf')
cairosvg.svg2png(url='$f', write_to='figs/png600/$b.png', dpi=600, scale=600/96)"
  cp -f "$f" "figs/svg/$b.svg"
  printf "  %-46s eps %-6s svg (source)\n" "$b" \
         "$(stat -c %s "figs/eps/$b.eps" | numfmt --to=iec)"
done

echo
echo "eps: $(ls -1 figs/eps/*.eps | wc -l)   svg: $(ls -1 figs/svg/*.svg | wc -l)   png600: $(ls -1 figs/png600/*.png | wc -l)"
echo "Not covered: Figure 2 (the algorithm pseudocode), which is not script-drawn."
