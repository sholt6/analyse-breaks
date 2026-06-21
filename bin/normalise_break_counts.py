#!/usr/bin/env python3

import argparse
import logging
import csv
from pathlib import Path


# Set up argument parser
def positive_int(value):
    value = int(value)
    if value <= 0:
        raise argparse.ArgumentTypeError(
            f"{value} is not a valid positive integer"
        )
    return value


parser = argparse.ArgumentParser(
    prog="normalise_break_counts.py",
    description="""Accepts a BED file containing break sites identified by
    INDUCE-Seq along with metadata and the same BED after filtering for 
    canonical restriction enzyme cutting sites. Outputs a normalised count of
    canonical breaks.""",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

parser.add_argument('raw_bed',
                    help="""path to an INDUCE-seq breaks BED""")

parser.add_argument('canonical_bed',
                    help="""path to BED filtered for canonical cut sites""")

parser.add_argument('-o', '--output-file', default="sample.counts.tsv",
                    help="""desired path to output TSV""")


parser.add_argument('-n', '--normalisation-value', type=positive_int, 
                    default=1000,
                    help="""raw sites count is divided by this positive 
                    integer value to give more readable results""")

parser.add_argument('-r', '--rounding-value', type=positive_int, 
                    default=3,
                    help="""normalised counts will be rounded to this many
                    decimal places""")

parser.add_argument('-l', '--log-file', default='mapq_filter_reads.log',
                    help="""path to log file""")


def get_sample_name(raw_bed: str) -> str:
    return raw_bed.split("/")[-1].split('.')[0]


def count_lines(file_name: str) -> int:
    return sum(1 for _ in open(file_name))


def normalise_counts(raw_sites_count: int, canonical_sites_count: int, 
                     normalisation_value: int, rounding_value: int) -> float:
    
    if raw_sites_count == 0:
        return 0
    
    return round( 
        ( canonical_sites_count / ( raw_sites_count / normalisation_value ) ), 
        rounding_value )


def write_results(sample_name: str, raw_sites_count: int,
                  canonical_sites_count: int, normalised_count: float,
                  output_file: str) -> None:

    with open(output_file, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        writer.writerow([sample_name, raw_sites_count,
                         canonical_sites_count, normalised_count])


# Functions
def main() -> None:
    args = parser.parse_args()

    logging.basicConfig(filename=args.log_file, level=logging.INFO, 
                        format="[%(asctime)s] %(levelname)s - %(message)s")
    logger = logging.getLogger(Path(__file__).stem)
    
    logger.info(f"Starting count normalisation with args: {args}")

    raw_bed = args.raw_bed
    canonical_bed = args.canonical_bed
    output_file = args.output_file
    normalisation_value = args.normalisation_value
    rounding_value = args.rounding_value

    sample_name = get_sample_name(raw_bed)

    raw_sites_count = count_lines(raw_bed)
    canonical_sites_count = count_lines(canonical_bed)

    if raw_sites_count == 0:
        logger.warning(f"No entries in raw breaks BED")

    normalised_count = normalise_counts(raw_sites_count, canonical_sites_count,
                                        normalisation_value, rounding_value)

    write_results(sample_name, raw_sites_count, canonical_sites_count, 
                  normalised_count, output_file)


    logger.info(f"Normalisation completed and written to {output_file}")

if __name__ == "__main__":
    main()
