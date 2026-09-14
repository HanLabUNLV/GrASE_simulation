"""Per-event provenance for merged bipartition counts.

PERMANENT PIPELINE STEP -- not a workaround. Required before
`infer_bipartition_gt.R --merged`, which will stop and name this script if the
output is missing.

merge_exon_sj_counts.R records, per bipartition side, whether the count used is
the exonic distinct set ("exon") or the side's exclusive junction ("sj"), plus
the junction ids. Both are constant per (gene,event). They are stored in two
different places ON PURPOSE:

  diff1_source / diff2_source   ->  carried per ROW in exontest.R's annotated
      output. Four characters, and add_significant() needs the source on every
      row to choose between --min_dpi and --min_dpi_sj.

  intron_distinct1 / intron_distinct2  ->  HERE, one row per (gene,event).
      These average 123 characters. The annotated output repeats each event
      about 14 times on a 10-contrast run, so carrying them per row would add
      1.47 GB to DICE's 4.30 GB TSS/TTS file to store 422,081 facts fourteen
      times. infer_bipartition_gt.R is their only consumer and needs them once
      per event, which is the shape this script writes.

Do not "simplify" by folding the ids into the annotated output; that trade was
measured and rejected. The same repetition is why a large merged count set is
trimmed of intron_distinct* before combining, which means the ids may exist
only in the per-gene files this script reads.

Usage: python3 scripts/extract_merge_provenance.py <countdir> <src> <out.txt>
  e.g. python3 scripts/extract_merge_provenance.py bipartition.merged.counts internal \\
         results/gt/merge_provenance.internal.txt
Out:  src gene event diff1_source diff2_source intron_distinct1 intron_distinct2
"""
import sys, glob, os, csv

countdir, src, out = sys.argv[1], sys.argv[2], sys.argv[3]
COLS = ["diff1_source", "diff2_source", "intron_distinct1", "intron_distinct2"]

n = 0
with open(out, "w") as fh:
    fh.write("\t".join(["src", "gene", "event"] + COLS) + "\n")
    for fp in sorted(glob.glob(os.path.join(countdir, "*.bipartition.exoncnt.txt"))):
        seen = set()
        with open(fp) as g:
            for r in csv.DictReader(g, delimiter="\t"):
                k = (r["gene"], r["event"])
                if k in seen:
                    continue
                seen.add(k)
                fh.write("\t".join([src, r["gene"], r["event"]] +
                                   [(r.get(c) or "NA") for c in COLS]) + "\n")
                n += 1
print(f"wrote {n} events -> {out}")
