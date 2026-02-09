#!/bin/bash
set -eux


CELLTYPE=$1
# Set input and output directories

# Build high-confidence SJ table:
# - canonical motifs only
# - appears in >=5 samples
# - >=10 unique reads total
# - >=10bp max overhang
awk '
{
  if ($5 != 1 && $5 != 2) next   # canonical only

  key = $1 FS $2 FS $3 FS $4
  samples[key]++
  uniq[key] += $7
  motif[key] = $5
  annot[key] = $6
  maxOH[key] = (maxOH[key] > $9 ? maxOH[key] : $9)
}
END {
  for (k in samples)
    if (samples[k] >= 5 && uniq[k] >= 10 && maxOH[k] >= 10) {
      split(k,a,FS)
      print a[1],a[2],a[3],a[4],motif[k],annot[k],uniq[k],maxOH[k]
    }
}' OFS="\t" ${CELLTYPE}.raw.SJ.tab > ${CELLTYPE}.highconf.SJ.tab

awk '
$7 >= 20 &&      # >=20 uniquely mapped reads
$8 >= 10          # >=10bp overhang
' ${CELLTYPE}.highconf.SJ.tab > ${CELLTYPE}.final.SJ.tab

# Count junctions
wc -l ${CELLTYPE}.final.SJ.tab

# Motif distribution
cut -f5 ${CELLTYPE}.final.SJ.tab | sort | uniq -c

# Annotated vs novel
cut -f6 ${CELLTYPE}.final.SJ.tab | sort | uniq -c

