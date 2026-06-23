#!/usr/bin/env python3

import argparse
import logging
import time
import csv
from pathlib import Path
from typing import Callable
from dataclasses import dataclass


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

parser.add_argument('-l', '--log-file', default='normalise_break_counts.log',
                    help="""path to log file""")


@dataclass
class NormalisationResult:
    sample_name: str=""
    raw_sites_count: int=0
    canonical_sites_count: int=0
    normalised_count: float=0.0

    def summary(self):
        return [self.sample_name, 
                self.raw_sites_count, 
                self.canonical_sites_count,
                self.normalised_count]


# Functions
def time_log(logger_name: str) -> Callable:
    """Decorator which measures duration of a function call and adds directly to log file"""
    def wrapper(func: Callable) -> Callable:
        def wrapper_func(*args, **kwargs):
            logger = logging.getLogger(logger_name)
            start = time.perf_counter()
            result = func(*args, **kwargs)
            end = time.perf_counter()
            logger.info(f"Function {func.__name__} completed in {end - start:.4f} seconds")
            return result
        return wrapper_func
    return wrapper


@time_log("logger")
def get_sample_name(raw_bed: str) -> str:
    """Extracts sample name from file name, expected format: <sample-name>.breakends.bed"""
    return raw_bed.split("/")[-1].split('.')[0]


@time_log("logger")
def count_lines(file_name: str) -> int:
    """Counts lines in given file and enables timing of this operation"""
    return sum(1 for _ in open(file_name))


@time_log("logger")
def normalise_counts(raw_sites_count: int, canonical_sites_count: int, 
                     normalisation_value: int, rounding_value: int) -> float:
    """Normalises count of canonical sites against total sites and renders more readable"""    
    if raw_sites_count == 0:
        return 0
    
    return round( 
        ( canonical_sites_count / ( raw_sites_count / normalisation_value ) ), 
        rounding_value )


@time_log("logger")
def write_results(result: NormalisationResult, output_file: str) -> None:
    """Writes results to specified output file"""
    with open(output_file, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        writer.writerow([result.sample_name, 
                         result.raw_sites_count,
                         result.canonical_sites_count, 
                         result.normalised_count])


@time_log("logger")
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

    result = NormalisationResult()

    result.sample_name = get_sample_name(raw_bed)

    result.raw_sites_count = count_lines(raw_bed)
    result.canonical_sites_count = count_lines(canonical_bed)

    if result.raw_sites_count == 0:
        logger.warning(f"No entries in raw breaks BED: {raw_bed}")

    result.normalised_count = normalise_counts(result.raw_sites_count, 
                                               result.canonical_sites_count,
                                               normalisation_value, 
                                               rounding_value)

    write_results(result, output_file)

    logger.info(f"Normalisation completed and written to {output_file}")


if __name__ == "__main__":
    main()
