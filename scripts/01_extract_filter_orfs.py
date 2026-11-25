#!/usr/bin/env python3
"""
trim_to_cds.py

Trims each HCV sequence in a fasta file to its polyprotein coding region,
by finding the longest open reading frame (ATG...stop) in each sequence
individually. This avoids relying on fixed coordinates (e.g. H77 numbering),
which don't apply uniformly across a mixed query+reference fasta file.

Sequences where no clean long ORF is found are written to a separate
"flagged" file for manual inspection, rather than silently mistrimmed.

Requirements: biopython  (pip install biopython)

Usage:
    python trim_to_cds.py input.fasta output_trimmed.fasta [--min-orf 8000]

Notes:
- HCV's polyprotein ORF is ~9024-9099 nt depending on genotype, so
  --min-orf defaults to 8000 as a safe lower bound that still excludes
  spurious short ORFs, but you can lower it if you're working with
  partial genome sequences rather than full-length genomes.
- This finds the ORF in the forward strand only, frames 0/1/2. HCV is a
  positive-sense ssRNA virus, so reverse-complement ORFs aren't expected;
  if any turn up, they're a sign that sequence's submitted orientation
  should be checked.
"""

import argparse

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
except ImportError:
    class SeqRecord:
        def __init__(self, seq, id="", description=""):
            self.seq = seq
            self.id = id
            self.description = description

        def __getitem__(self, index):
            return SeqRecord(self.seq[index], id=self.id, description=self.description)

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

        @staticmethod
        def write(records, handle_or_path, fmt):
            if fmt != "fasta":
                raise ValueError(f"Unsupported format: {fmt}")
            is_path = isinstance(handle_or_path, str)
            f = open(handle_or_path, "w", encoding="utf-8") if is_path else handle_or_path
            try:
                for r in records:
                    f.write(f">{r.description}\n")
                    seq_str = str(r.seq)
                    for i in range(0, len(seq_str), 80):
                        f.write(seq_str[i:i+80] + "\n")
            finally:
                if is_path:
                    f.close()

STOP_CODONS = {"TAA", "TAG", "TGA"}


def find_longest_orf(seq_str, min_len):
    """Return (start, end) 0-based, end exclusive, of the longest ATG...stop
    ORF across the 3 forward reading frames. Returns None if nothing found
    at least min_len long."""
    seq_str = seq_str.upper().replace("U", "T")
    best = None
    for frame in range(3):
        i = frame
        while i < len(seq_str) - 2:
            codon = seq_str[i:i + 3]
            if codon == "ATG":
                # scan forward for the first in-frame stop
                j = i
                while j < len(seq_str) - 2:
                    c = seq_str[j:j + 3]
                    if c in STOP_CODONS:
                        orf_len = j + 3 - i
                        if orf_len >= min_len and (best is None or orf_len > best[1] - best[0]):
                            best = (i, j + 3)
                        break
                    j += 3
                else:
                    # ran off the end without a stop -- sequence may be
                    # incomplete at the 3' end; record as candidate anyway
                    orf_len = len(seq_str) - i
                    if orf_len >= min_len and (best is None or orf_len > best[1] - best[0]):
                        best = (i, len(seq_str))
            i += 3
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input_fasta")
    ap.add_argument("output_fasta")
    ap.add_argument("--min-orf", type=int, default=8000,
                     help="minimum ORF length in nt to accept (default 8000; "
                          "lower this for partial-genome sequences)")
    ap.add_argument("--trim-incomplete-stop", action="store_true",
                   help="trim a trailing partial stop codon from an ORF reaching the sequence end")
    args = ap.parse_args()

    kept, flagged = [], []

    for rec in SeqIO.parse(args.input_fasta, "fasta"):
        seq_str = str(rec.seq)
        orf = find_longest_orf(seq_str, args.min_orf)
        if orf is None:
            flagged.append(rec)
            continue
        start, end = orf
        remainder = (end - start) % 3
        if args.trim_incomplete_stop and end == len(seq_str) and remainder:
            tail = seq_str[end - remainder:end].upper().replace("U", "T")
            if any(stop.startswith(tail) for stop in STOP_CODONS):
                end -= remainder
        trimmed = rec[start:end]
        trimmed.id = rec.id
        trimmed.description = f"{rec.description} | trimmed {start}:{end} (len={end-start})"
        kept.append(trimmed)

    SeqIO.write(kept, args.output_fasta, "fasta")
    flagged_path = args.output_fasta.rsplit(".", 1)[0] + "_FLAGGED.fasta"
    if flagged:
        SeqIO.write(flagged, flagged_path, "fasta")

    print(f"Trimmed OK:  {len(kept)} sequences -> {args.output_fasta}")
    print(f"Flagged:     {len(flagged)} sequences -> {flagged_path if flagged else '(none)'}")
    if flagged:
        print("\nFlagged sequences had no ORF >= --min-orf nt in any forward frame.")
        print("Common causes: partial/fragment sequence, non-standard start,")
        print("sequencing gaps (N runs) breaking the ORF, or wrong orientation.")
        print("Inspect these manually -- do not assume they're unusable.")


if __name__ == "__main__":
    main()
