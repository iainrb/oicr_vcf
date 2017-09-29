#! /usr/bin/env python

# script to parse VCF files, extract basic stats and write in JSON format

# separate JSON file for each individual

import re, sys, json

def main():

    stats = parse_stats(sys.stdin) # data structure
    # stats has results for all individuals; write as JSON for each individual


def parse_stats(infile):

    # want to read header (as a whole) and body (in chunks)
    # need stats for each individual
    stats = {}

    buffer_size = 50 * 10**6 # input buffer size, in bytes

    input_lines = infile.readlines(buffer_size)
    
    # read VCF header

    input_lines = infile.readlines(buffer_size)
    header_lines = []
    column_heads = None
    for line in input_lines:
        if re.match('##', line):
            header_lines.append(line)
        elif re.match('#CHROM', line):
            column_heads = line
            break
        else:
            msg = "Unexpected line in VCF header; line does not start "+\
                  "with '##' or '#CHROM': "+line
            raise ValueError(msg)
    # next line will be first line in VCF body

    # TODO warn when input exceeds eg. 1 TB?
    
    # read VCF body
    while True:
        lines = infile.readlines(buffer_size)
        if lines == []: break
        for line in lines:
            fields = parse_body_line(line)
            stats = update_stats(fields)

    return stats
            

if __name__ == "__main__":
    main()
