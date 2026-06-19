#!/usr/bin/env python

import sys
import csv

def bed_reader(file_name):
    with open(file_name, "r") as file:
        reader = csv.reader(file, delimiter="\t")
        for row in reader:
            yield row


def filter_by_score(rows, min_score=30):
    for row in rows:
        if int(row[4]) >= 30:
            yield row


def main():
    input_bed = sys.argv[1]
    output_bed = sys.argv[2]
    mapq_threshold = sys.argv[3]

    read_bed = bed_reader(input_bed)
    filtered_bed_rows = filter_by_score(read_bed, mapq_threshold)
    
    with open('filtered.bed', "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        for row in filtered_bed_rows:
            writer.writerow(row)

if __name__ == "__main__":
    main()
