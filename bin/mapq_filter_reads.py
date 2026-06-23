#!/usr/bin/env python3

import argparse
import logging
import time
import csv
from pathlib import Path
from typing import Callable, Generator
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
def bed_reader(file_name: str) -> Generator[list[str], None, None]:
    """Generator for more memory-safe BED file reading"""
    with open(file_name, "r") as file:
        reader = csv.reader(file, delimiter="\t")
        for row in reader:
            yield row


@time_log("logger")
def filter_bed(read_bed: Generator[list[str], None, None], mapq_threshold: int,
                filtered_bed: str) -> None:
    """Iterates over bed_reader and filters on specified MapQ.
    Outputs directly to filtered BED"""
    with open(filtered_bed, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        for row in (row for row in read_bed if int(row[4]) >= mapq_threshold):
            writer.writerow(row)


@time_log("logger")
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

    filter_bed(read_bed, mapq_threshold, filtered_bed)

    logger.info(f"Filtering completed and written to {filtered_bed}")


if __name__ == "__main__":
    main()
