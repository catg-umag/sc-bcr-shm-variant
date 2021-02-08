#!/usr/bin/env python3
import argparse
import csv


def main():
    args = parse_arguments()

    with open(args.input) as inf, open(args.output, "w") as outf:
        reader = csv.DictReader(inf, delimiter="\t")
        writer = csv.writer(outf)

        # header
        writer.writerow(
            [
                "sequence_id",
                "vdj_start",
                "vdj_end",
                "cdr1_start",
                "cdr1_end",
                "cdr2_start",
                "cdr2_end",
                "cdr3_start",
                "cdr3_end",
            ]
        )

        columns = [
            "sequence_id",
            "v_sequence_start",
            "j_sequence_end",
            "cdr1_start",
            "cdr1_end",
            "cdr2_start",
            "cdr2_end",
            "cdr3_start",
            "cdr3_end",
        ]
        for row in reader:
            writer.writerow([row[x] for x in columns])


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Extract VDJ and CDR regions from Change-O Igblast AIRR output"
    )
    parser.add_argument(
        "--input", "-i", required=True, help="Change-O AIRR output (tsv)"
    )
    parser.add_argument("--output", "-o", required=True, help="output file (csv)")

    return parser.parse_args()


if __name__ == "__main__":
    main()