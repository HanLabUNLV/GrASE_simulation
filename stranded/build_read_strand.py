"""Recover each simulated read's TRUE strand of origin.

simulate_reads.R calls polyester's simulate_experiment() without
strand_specific, so the library is unstranded: mate orientation is random and
carries no information about the source gene's strand (checked directly -- mate1
is ~60/40 forward in both a + strand and a - strand gene body). dexseq_count.sh
nevertheless ran with --stranded reverse, so the existing counts were filtered
by an orientation rule that means nothing: it drops about half the reads at
random and does NOT stop an antisense neighbour's reads landing in a gene's bins.

polyester does write the source transcript into the read name

    >read33488787/ENST00000361624.2;mate1:402-501;mate2:510-609

and STAR truncates it at the first "/", so the BAM keeps only "read33488787".
The mapping survives in the FASTA. Transcript -> gene -> strand then gives the
read its true strand, which is exactly what a stranded protocol would have
delivered. Note this recovers STRAND ONLY, not gene: same-strand overlaps stay
ambiguous, as they would in a real stranded experiment.

Read ids are "read<N>" with N < ~4e7, so a flat int8 array indexed by N is the
cheapest lookup (0 = unknown, 1 = +, 2 = -).

Usage: python3 scripts/build_read_strand.py <mate1.fa.gz> <out.npy>
"""
import sys, gzip, re, numpy as np

BASE = "/mnt/data1/home/mirahan/GrASE_simulation"
GTF  = f"{BASE}/ref/gencode.v28.annotation.gtf"

def tx_strand():
    """transcript_id -> '+'/'-' from the reference GTF"""
    out = {}
    pat = re.compile(r'transcript_id "([^"]+)"')
    with open(GTF) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 9)
            if len(f) < 9 or f[2] != "transcript":
                continue
            m = pat.search(f[8])
            if m:
                out[m.group(1)] = f[6]
    return out

def main(fa, out):
    ts = tx_strand()
    sys.stderr.write(f"transcripts in GTF: {len(ts)}\n")
    arr = np.zeros(50_000_000, dtype=np.int8)
    n = miss = maxid = 0
    op = gzip.open if fa.endswith(".gz") else open
    with op(fa, "rt") as fh:
        for line in fh:
            if line[0] != ">":
                continue
            slash = line.index("/")
            rid = int(line[5:slash])            # skip ">read"
            semi = line.index(";", slash)
            st = ts.get(line[slash + 1:semi])
            if st is None:
                miss += 1
            else:
                arr[rid] = 1 if st == "+" else 2
            n += 1
            if rid > maxid:
                maxid = rid
    sys.stderr.write(f"reads {n}, max id {maxid}, transcript not in GTF {miss}\n")
    np.save(out, arr[:maxid + 1])
    sys.stderr.write(f"wrote {out}\n")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
