#!/usr/bin/env python3
import argparse
from typing import Dict, List
import pandas as pd


def main():
    args = parse_arguments()

    df = (
        pd.read_csv(args.input)
        .query("ref_vdj_coverage == 1")
        .drop(columns=["consensus"])
    )

    selected_cells = (
        df.groupby("cell")
        .size()
        .to_frame("n")
        .query(f"n >= {args.min_umis}")
        .index.values
    )

    variants = []
    for cell in selected_cells:
        subdf = df[df.cell == cell]
        umis = subdf.umi
        sequences = subdf.aligned_consensus.values
        depths = subdf.depths.values

        for i in range(len(sequences[0])):
            cons = set(s[i] for s in sequences) - set(("-", "N"))
            if len(cons) > 1:
                counts = count_bases([s[i] for s in sequences])
                ratios = [
                    x / sum(counts.values())
                    for x in sorted(counts.values(), reverse=True)
                ]

                cons_nucl = "N"
                for n, c in counts.items():
                    if c == ratios[0] * sum(counts.values()):
                        cons_nucl = n
                        break

                if (
                    sum(counts.values()) >= args.min_umis
                    and ratios[1] >= args.min_ratio
                ):
                    for (umi, seq) in zip(umis, sequences):
                        if seq[i] not in ("N", "-"):
                            variants.append(
                                {
                                    "cell": cell,
                                    "umi": umi,
                                    "pos": i + 1,
                                    "cons_nucl": cons_nucl,
                                    "nucl": seq[i],
                                    "depth": depths[i]
                                }
                            )

    pd.DataFrame(variants).to_csv(args.output, index=False)


def count_bases(bases: List[str]) -> Dict[str, int]:
    d = {"A": 0, "C": 0, "G": 0, "T": 0}
    for b in bases:
        if b in d.keys():
            d[b] += 1

    return d


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Search places of possible SHM occurence"
    )
    parser.add_argument(
        "--input", "-i", required=True, help="Consensus list by Cell/UMI (.csv)"
    )
    parser.add_argument("--output", "-o", required=True, help="output file (.csv)")
    parser.add_argument(
        "--min-umis",
        "-u",
        default=2,
        help="Minimun ammount of UMIs supporting cell",
        dest="min_umis",
    )
    parser.add_argument(
        "--min-ratio",
        "-r",
        default=0,
        help="Minimun alt allele ratio to consider it a valid case",
        dest="min_ratio",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()