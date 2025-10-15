#!/data/stahlberg/Manuel/envs/umierrorGYN/bin/python
'''
UMI error correct, split_bam.py - Separate aligned reads based on genomic coordiantes into panels A and B.

:Author: Manuel Luna

Purpose
-------

Separate aligned reads based on genomic coordiantes into panels A and B.

'''
import pysam
import os
import argparse
import logging
from collections import defaultdict
from datetime import datetime

def parse_args():
    p = argparse.ArgumentParser(
        description="Split BAM reads into A/B/dump using two BEDs of amplicons and start-pad rules."
    )
    p.add_argument("-b", "--bam", dest="bam_file", required=True,
                   help="Input BAM (aligned; index not required for streaming).")
    p.add_argument("--bed_a", required=True,
                   help="BED for group A: chrom start end name (0-based, end-exclusive).")
    p.add_argument("--bed_b", required=True,
                   help="BED for group B: chrom start end name (0-based, end-exclusive).")
    p.add_argument("-o", "--output_dir", required=True, help="Output directory.")
    p.add_argument("-s", "--sample", required=True, help="Sample name/prefix for outputs.")
    p.add_argument("--min_start_pad", type=int, default=12,
                   help="Min nt read start before amplicon start (default 12).")
    p.add_argument("--max_start_pad", type=int, default=32,
                   help="Max nt read start before amplicon start (default 32).")
    p.add_argument("--max_end_pad", type=int, default=100,
                   help="Max nt the read may extend beyond amplicon end (default 100).")
    return p.parse_args()

def read_bed(path, label):
    """Return chrom->list of {'name','start','end','label'}."""
    idx = defaultdict(list)
    with open(path) as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end   = int(parts[2])
            name  = parts[3] if len(parts) > 3 else f"{chrom}:{start}-{end}"
            idx[chrom].append({"name": name, "start": start, "end": end, "label": label})
    return idx

def read_ref_positions_including_deletions(read):
    """
    Return a set of reference positions covered by the alignment.
    Includes matches/mismatches and deletions (qpos=None,rpos!=None).
    """
    ref_pos = set()
    for qpos, rpos in read.get_aligned_pairs(matches_only=False):
        if rpos is not None:
            ref_pos.add(rpos)
    return ref_pos
    
def choose_amplicon_for_read(read, bedA, bedB, min_start_pad, max_start_pad, max_end_pad):
    """
    Decide A/B assignment for a mapped read:
      - Must fully cover the amplicon [start,end) (deletions count as coverage).
      - start_delta = a.start - read.reference_start must be within [min_start_pad, max_start_pad].
      - (Optional) end_delta = (read.reference_end-1) - (a.end-1) <= max_end_pad.
      - Prefer inside-window candidates; then smaller start_delta; then shorter amplicon; then A over B.
    Return dict {'name','start','end','label'} or None.
    """
    chrom = read.reference_name
    if chrom is None:
        return None

    listA = bedA.get(chrom, [])
    listB = bedB.get(chrom, [])
    if not listA and not listB:
        return None

    rstart = read.reference_start
    rend   = read.reference_end  # exclusive
    ref_set = read_ref_positions_including_deletions(read)
    if not ref_set:
        return None

    candidates = []
    for amp in listA + listB:
        a0, a1 = amp["start"], amp["end"]  # [a0,a1)
        # Must fully cover amplicon:
        if not all((pos in ref_set) for pos in range(a0, a1)):
            continue

        start_delta = a0 - rstart
        in_window = (min_start_pad <= start_delta <= max_start_pad)
        end_delta = (rend - 1) - (a1 - 1)
        if end_delta > max_end_pad:
            continue

        if in_window:
            score = (0, start_delta, (a1 - a0))  # inside-window first; prefer smaller delta; then shorter amplicon
        else:
            # distance to window if outside
            if start_delta < min_start_pad:
                dist_to_window = min_start_pad - start_delta
            else:
                dist_to_window = start_delta - max_start_pad
            score = (1, dist_to_window, (a1 - a0))

        candidates.append((score, amp))

    if not candidates:
        return None

    # Sort; final tie-breaker prefers A over B
    candidates.sort(key=lambda x: (x[0][0], x[0][1], x[0][2], 0 if x[1]["label"] == "A" else 1))
    return candidates[0][1]

def split_bam(output_dir,sample,bamfile,bed_a,bed_b, min_start_pad=12, max_start_pad=32, max_end_pad=100):
    outA_unsorted = os.path.join(output_dir, f"{sample}_A.unsorted.bam")
    outB_unsorted = os.path.join(output_dir, f"{sample}_B.unsorted.bam")
    outD_unsorted = os.path.join(output_dir, f"{sample}_dump.unsorted.bam")

    outA = os.path.join(output_dir, f"{sample}_A.sorted.bam")
    outB = os.path.join(output_dir, f"{sample}_B.sorted.bam")
    outD = os.path.join(output_dir, f"{sample}_dump.sorted.bam")

    log_summary = os.path.join(output_dir, f"{sample}_split_summary.tsv")

    logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(message)s')
    logging.info("Started BAM splitting by two BEDs")
    logging.info(f"Sample: {sample} | BAM: {bamfile}")

    # Load BEDs
    bedA = read_bed(bed_a, label="A")
    bedB = read_bed(bed_b, label="B")
    nA = sum(len(v) for v in bedA.values())
    nB = sum(len(v) for v in bedB.values())
    logging.info(f"Loaded {nA} amplicons from A and {nB} from B")

    counts = {"A": 0, "B": 0, "dump": 0}

    # Stream input and route reads
    with pysam.AlignmentFile(bamfile, "rb") as infile, \
         pysam.AlignmentFile(outA_unsorted, "wb", template=infile) as outA_fh, \
         pysam.AlignmentFile(outB_unsorted, "wb", template=infile) as outB_fh, \
         pysam.AlignmentFile(outD_unsorted, "wb", template=infile) as outD_fh:

        for read in infile:
            if read.is_unmapped:
                outD_fh.write(read)
                counts["dump"] += 1
                continue

            amp = choose_amplicon_for_read(
                read, bedA, bedB,
                min_start_pad, max_start_pad, max_end_pad
            )

            if amp is None:
                outD_fh.write(read)
                counts["dump"] += 1
            else:
                if amp["label"] == "A":
                    outA_fh.write(read)
                    counts["A"] += 1
                else:
                    outB_fh.write(read)
                    counts["B"] += 1

    # Sort & index outputs
    logging.info("Sorting and indexing outputs...")
    for src, dst in [(outA_unsorted, outA), (outB_unsorted, outB), (outD_unsorted, outD)]:
        pysam.sort("-o", dst, src)
        pysam.index(dst)
        os.remove(src)

    # Write summary
    logging.info("Writing summary...")
    with open(log_summary, "w") as f:
        f.write("Bin\tReads\n")
        f.write(f"A\t{counts['A']}\n")
        f.write(f"B\t{counts['B']}\n")
        f.write(f"dump\t{counts['dump']}\n")
    return outA, outB
    logging.info("Done.")

if __name__ == "__main__":
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    split_bam(args.output_dir,args.sample,args.bamfile,args.bed_a,args.bed_b,args.min_start_pad, args.max_start_pad, args.max_end_pad)
