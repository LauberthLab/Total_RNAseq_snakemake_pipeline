#!/usr/bin/env python
"""Tidy a featureCounts output file.

- Drops Chr/Start/End/Strand columns.
- Renames sample columns from BAM paths to clean sample ids
  (handles trailing "_S<n>" lane tags and STAR-style suffixes).
"""

import os
import sys
import pandas as pd


def extract_sample_name(bam_path: str) -> str:
    base = os.path.basename(bam_path)
    if "Aligned.sortedByCoord.out.bam" in base:
        return base.replace("Aligned.sortedByCoord.out.bam", "")
    base = base.replace(".bam", "")
    return base


def main() -> None:
    if len(sys.argv) != 3:
        print("Usage: process_fc.py <input_file> <output_file>", file=sys.stderr)
        sys.exit(1)

    input_file, output_file = sys.argv[1], sys.argv[2]

    fc = pd.read_csv(input_file, sep="\t", header=1)
    fc = fc.drop(columns=["Chr", "Start", "End", "Strand"], errors="ignore")
    fc.columns = [
        c if c in ("Geneid", "Length") else extract_sample_name(c)
        for c in fc.columns
    ]

    out_dir = os.path.dirname(os.path.abspath(output_file))
    os.makedirs(out_dir, exist_ok=True)
    fc.to_csv(output_file, sep="\t", index=False)

    print(f"Processed {len(fc)} features x {len(fc.columns) - 2} samples -> {output_file}")


if __name__ == "__main__":
    main()