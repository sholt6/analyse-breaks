#!/usr/bin/env python3

import csv
import argparse
import logging
from pathlib import Path
from collections.abc import Generator


# Set up argument parser
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

parser.add_argument('-q', '--mapq-threshold', type=int, default=30,
                    help="""path to an INDUCE-seq breaks BED""")

parser.add_argument('-l', '--log-file', default='mapq_filter_reads.log',
                     help="""path to log file""")


# Functions
def bed_reader(file_name: str) -> Generator[list[str], None, None]:
    with open(file_name, "r") as file:
        reader = csv.reader(file, delimiter="\t")
        for row in reader:
            yield row


def main() -> None:
    args = parser.parse_args()

    logging.basicConfig(filename=args.log_file, level=logging.INFO, 
                        format="[%(asctime)s] %(levelname)s - %(message)s")
    logger = logging.getLogger(Path(__file__).stem)

    logger.info(f"Starting filtering with args: {args}")

    input_bed = args.input_bed
    filtered_bed = args.filtered_bed
    mapq_threshold = args.mapq_threshold

    read_bed = bed_reader(input_bed)

    filtered_bed_rows = ( row for row in read_bed if int(row[4]) >= mapq_threshold )

    with open(filtered_bed, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        for row in filtered_bed_rows:
            writer.writerow(row)

    logger.info(f"Filtering completed and written to {filtered_bed}")

if __name__ == "__main__":
    main()
