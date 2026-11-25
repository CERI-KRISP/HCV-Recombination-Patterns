#!/usr/bin/env python3
"""Map GARD partition boundaries using an archived retained-column TSV."""

import argparse
import csv
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("coordinate_map", type=Path)
    parser.add_argument("left_columns", type=int, nargs="+",
                        help="One-based GARD alignment columns to the left of each boundary.")
    args = parser.parse_args()

    with args.coordinate_map.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        source = {
            int(row["gard_alignment_column_1based"]): int(row["source_alignment_column_1based"])
            for row in reader
        }

    if any(n < 1 or n not in source or n + 1 not in source for n in args.left_columns):
        parser.error("each boundary must have both flanks in the coordinate map")

    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(("gard_left_column", "source_left_column", "source_right_column"))
    for n in args.left_columns:
        writer.writerow((n, source[n], source[n + 1]))


if __name__ == "__main__":
    main()
