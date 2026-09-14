"""Split an alignment BAM by each read's TRUE strand of origin.

The reads are unstranded (polyester default), so their FLAG orientation says
nothing about the source gene's strand. scripts/build_read_strand.py recovers
the real strand from the polyester read name. This script carries that label
into the alignments by writing two BAMs -- reads from + strand transcripts and
reads from - strand transcripts.

Counted against a strand-split GFF with --stranded no, a plus-origin read can
then only reach a plus-strand gene. That is exactly the discrimination a
stranded protocol provides, without rewriting any FLAG (which would leave SEQ
and CIGAR inconsistent with the claimed strand). Strand only, never gene:
same-strand overlaps stay ambiguous, as they would in a real experiment.

Usage: python3 scripts/split_bam_by_origin.py <in.bam> <strand.npy> <out_prefix>
Out:   <out_prefix>.plus.bam(.bai), <out_prefix>.minus.bam(.bai)
"""
import sys, numpy as np, pysam

def main(bam_in, npy, prefix):
    arr = np.load(npy)
    n_ids = len(arr)
    src = pysam.AlignmentFile(bam_in, "rb")
    op = pysam.AlignmentFile(prefix + ".plus.bam",  "wb", template=src)
    om = pysam.AlignmentFile(prefix + ".minus.bam", "wb", template=src)
    n = np_ = nm = unk = 0
    for r in src.fetch(until_eof=True):
        n += 1
        try:
            rid = int(r.query_name[4:])          # "read33488787" -> 33488787
        except ValueError:
            unk += 1; continue
        s = arr[rid] if rid < n_ids else 0
        if s == 1:
            op.write(r); np_ += 1
        elif s == 2:
            om.write(r); nm += 1
        else:
            unk += 1
    src.close(); op.close(); om.close()
    pysam.index(prefix + ".plus.bam")
    pysam.index(prefix + ".minus.bam")
    sys.stderr.write(f"{bam_in}: {n} records -> plus {np_}, minus {nm}, unknown {unk}\n")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
