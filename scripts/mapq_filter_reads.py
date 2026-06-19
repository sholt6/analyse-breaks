#!/usr/bin/env python

import csv
import argparse

parser = argparse.ArgumentParser(
    prog="mapq_filter_reads.py",
    description="""Accepts a BED file containing break sites identified by
    INDUCE-Seq along with metadata. Filters on quality column and outputs to a
    new BED file.""",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

parser.add_argument('input_bed',
                    help="""path to an INDUCE-seq breaks BED""")

parser.add_argument('filtered_bed',
                    help="""desired path to output filtered BED""")

parser.add_argument('mapq_threshold', type=int, default=30,
                    help="""path to an INDUCE-seq breaks BED""")


def bed_reader(file_name):
    with open(file_name, "r") as file:
        reader = csv.reader(file, delimiter="\t")
        for row in reader:
            yield row


def filter_by_score(rows, min_score=30):
    for row in rows:
        if int(row[4]) >= min_score:
            yield row


def main():
    args = parser.parse_args()

    input_bed = args.input_bed
    filtered_bed = args.filtered_bed
    mapq_threshold = args.mapq_threshold

    read_bed = bed_reader(input_bed)
    filtered_bed_rows = filter_by_score(read_bed, mapq_threshold)
    
    with open(filtered_bed, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        for row in filtered_bed_rows:
            writer.writerow(row)


if __name__ == "__main__":
    main()
