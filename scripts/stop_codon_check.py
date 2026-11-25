#!/usr/bin/env python3
"""
stop_codon_check.py

QC step for the HCV codon-aware alignment workflow. Scans the unaligned
amino acid translations produced by TranslatorX (<prefix>.aaseqs.fasta)
and reports, per sequence:
  - whether any INTERNAL stop codons (*) are present (a genuine problem --
    indicates a frameshift or a true premature stop in that record)
  - how many ambiguous-residue calls (X) are present, and at what
    percentage of sequence length (informational -- these come from IUPAC
    ambiguity codes in the nucleotide sequence, e.g. mixed-base Sanger
    calls, and are NOT stop codons)
  - a PASS / REVIEW / FAIL flag based on the thresholds below

Output: a single TSV, one row per sequence, suitable for archiving
alongside the manuscript/methods as a reproducible QC record.

Usage:
    python3 stop_codon_check.py <prefix>.aaseqs.fasta output_qc.tsv [--x-warn-pct 2.0]
"""

import argparse
import csv

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError:
    class SeqRecord:
        def __init__(self, seq, id="", description=""):
            self.seq = seq
            self.id = id
            self.description = description

    class SeqIO:
        @staticmethod
        def parse(handle_or_path, fmt):
            if fmt != "fasta":
                raise ValueError(f"Unsupported format: {fmt}")
            is_path = isinstance(handle_or_path, str)
            f = open(handle_or_path, "r", encoding="utf-8") if is_path else handle_or_path
            try:
                cur_id = None
                cur_desc = None
                seq_chunks = []
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if cur_id is not None:
                            yield SeqRecord("".join(seq_chunks), id=cur_id, description=cur_desc)
                        cur_desc = line[1:]
                        cur_id = cur_desc.split()[0]
                        seq_chunks = []
                    else:
                        seq_chunks.append(line)
                if cur_id is not None:
                    yield SeqRecord("".join(seq_chunks), id=cur_id, description=cur_desc)
            finally:
                if is_path:
                    f.close()


def classify(internal_stops, x_pct, x_warn_pct):
    if internal_stops > 0:
        return "FAIL"
    if x_pct >= x_warn_pct:
        return "REVIEW"
    return "PASS"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("aaseqs_fasta", help="TranslatorX unaligned AA translation file (<prefix>.aaseqs.fasta)")
    ap.add_argument("output_tsv", help="Path to write the QC report TSV")
    ap.add_argument("--x-warn-pct", type=float, default=2.0,
                     help="Flag a sequence as REVIEW if ambiguous (X) residues "
                          "exceed this percentage of sequence length (default 2.0)")
    args = ap.parse_args()

    rows = []
    for rec in SeqIO.parse(args.aaseqs_fasta, "fasta"):
        seq = str(rec.seq)
        length = len(seq)
        # Internal stop = '*' anywhere except the final residue, which is
        # the expected terminal stop codon of the ORF.
        body = seq[:-1] if seq.endswith("*") else seq
        internal_stops = body.count("*")
        internal_stop_positions = [i + 1 for i, c in enumerate(body) if c == "*"]
        x_count = seq.count("X")
        x_pct = 100 * x_count / length if length else 0.0
        flag = classify(internal_stops, x_pct, args.x_warn_pct)

        rows.append({
            "sequence_id": rec.id,
            "aa_length": length,
            "internal_stop_count": internal_stops,
            "internal_stop_positions_aa": ";".join(map(str, internal_stop_positions)) if internal_stop_positions else "",
            "ambiguous_X_count": x_count,
            "ambiguous_X_pct": round(x_pct, 2),
            "flag": flag,
        })

    with open(args.output_tsv, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "sequence_id", "aa_length", "internal_stop_count",
                "internal_stop_positions_aa", "ambiguous_X_count",
                "ambiguous_X_pct", "flag",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    n_fail = sum(1 for r in rows if r["flag"] == "FAIL")
    n_review = sum(1 for r in rows if r["flag"] == "REVIEW")
    n_pass = sum(1 for r in rows if r["flag"] == "PASS")

    print(f"QC report written to {args.output_tsv}")
    print(f"  Total sequences: {len(rows)}")
    print(f"  PASS:   {n_pass}")
    print(f"  REVIEW: {n_review}  (ambiguous-X density >= {args.x_warn_pct}%, no stop codons)")
    print(f"  FAIL:   {n_fail}  (genuine internal stop codon detected)")
    if n_fail:
        print("\nFAIL sequences (real internal stop codons -- do not silently drop; "
              "mask the offending codon in the original nt fasta and re-run):")
        for r in rows:
            if r["flag"] == "FAIL":
                print(f"  {r['sequence_id']}  (position(s) aa: {r['internal_stop_positions_aa']})")


if __name__ == "__main__":
    main()
