#! /usr/bin/env python

# script to parse VCF files, extract basic stats and write in JSON format

# separate JSON file for each individual

import re, sys, json


def main():

    parser = vcf_parser()
    stats = parser.parse_stats(sys.stdin)


class vcf_parser:


    def __init__(self):
        self.buffer_size = 50 * 10**6 # input buffer size, in bytes


    def parse_stats(self, infile):

        # want to read header (as a whole) and body (in chunks)
        # need stats for each individual
        stats = {}

        # read VCF header
        header_lines = []
        column_heads = None

        while True:
            # read header lines one at a time
            line = infile.readline()
            if line == '':
                msg = "Reached end of file without finding end of VCF header"
                raise ValueError(msg)
            if re.match('##', line):
                header_lines.append(line)
            elif re.match('#CHROM', line):
                column_heads = line
                break
            else:
                msg = "Unexpected line in VCF header; line "+\
                    "does not start with '##' or '#CHROM': "+line
                raise ValueError(msg)

        sys.stderr.write("Read "+str(len(header_lines))+\
                             " lines in VCF header\n")
        if column_heads != None:
            sys.stderr.write("Read VCF column headers\n")

        # next line will be start of VCF body

        # TODO warn when input exceeds eg. 1 TB?
    
        # read VCF body
        line_count = 0
        while True:
            lines = infile.readlines(self.buffer_size)
            if lines == []: break
            for line in lines:
                line_count += 1
                #fields = parse_body_line(line)
                #stats = update_stats(fields)
        sys.stderr.write("Read "+str(line_count)+" lines in VCF body\n")
        return stats
            

# end of vcf_parser

if __name__ == "__main__":
    main()
