#!/bin/bash
# Re-test n_choose_2 on the SAME stranded counts with --min_reads=0, so the
# cross comparison can be scored with no read floor on any structure.
#
# phi is REUSED from the min_reads=10 run: the dispersion is estimated from the
# counts, which are identical, and reusing it keeps the only difference between
# the two runs the floor itself. Re-estimating would confound the comparison
# and cost another several hours.
#
# internal runs first (~7 min); TSSTTS is the long one (~6.5 h).
#
# Usage: bash scripts/run_nc2_nofloor_tests.sh
set -u
cd /mnt/data1/home/mirahan/GrASE_simulation
S=/mnt/data1/home/mirahan/GrASE/Rpkg/scripts
OUT=n_choose_2.stranded.test.EBapprox.nofloor
SRC=n_choose_2.stranded.test.EBapprox
mkdir -p "$OUT"

# exontest.R resolves --phi RELATIVE TO --outdir and reuses the file if it is
# already there (see the file.exists(outdir/phifile) check). So the phi files
# must sit INSIDE $OUT under a BARE name -- passing a directory-qualified path
# makes it look for $OUT/$SRC/phi... which does not exist, and it silently
# re-estimates and then dies writing the result. Symlink rather than copy:
# the TSS/TTS phi and its moderated companion are ~420 MB.
for f in "$SRC"/phi.nc2.*.txt; do
  [ -e "$f" ] || continue
  ln -sf "$(cd "$SRC" && pwd)/$(basename "$f")" "$OUT/$(basename "$f")"
done
echo "linked $(ls "$OUT"/phi.nc2.*.txt 2>/dev/null | wc -l) phi files into $OUT"

for a in internal TSSTTS; do
  echo "=== n_choose_2 $a min_reads=0  $(date +%H:%M:%S) ==="
  Rscript $S/exontest.R \
    --file=n_choose_2.${a}.exoncnt.combined.txt \
    --countdir=n_choose_2.${a}.stranded.counts/ \
    --outdir="$OUT" \
    --splittype=n_choose_2 --model=betabinom_EBapprox \
    --phi=phi.nc2.${a}.stranded.txt \
    --use_phi_loess --independent_filtering \
    --min_reads=0 --cond1=group1 --cond2=group2 \
    &> nc2.${a}.stranded.nofloor.log
  echo "  -> $(stat -c %s $OUT/test_n_choose_2.${a}_betabinom_EBapprox.txt 2>/dev/null || echo FAILED)"
done
echo "=== nc2 no-floor tests done $(date +%H:%M:%S) ==="
