#! /usr/bin/env python

# script to parse VCF files, extract basic stats and write in JSON format

# separate JSON file for each individual

import re, sys, json


def main():

    parser = vcf_parser()
    stats = parser.parse_stats(sys.stdin)


class vcf_parser:


    def __init__(self):
        self.buffer_size = 1 * 10**6 # input buffer size, in bytes
        self.total_fields = None

    def parse_stats(self, infile):

        # want to read header (as a whole) and body (in chunks)
        # need stats for each individual
        stats = {}

        # read VCF meta lines and header
        meta_lines = []
        header = None

        while True:
            # read header lines one at a time
            line = infile.readline()
            if line == '':
                msg = "Reached end of file without finding end of VCF header"
                raise ValueError(msg)
            if re.match('##', line):
                meta_lines.append(line)
            elif re.match('#CHROM', line):
                header = line
                break
            else:
                msg = "Unexpected line in VCF header; line "+\
                    "does not start with '##' or '#CHROM': "+line
                raise ValueError(msg)

        sys.stderr.write("Read "+str(len(meta_lines))+\
                             " lines in VCF metadata\n")
        (total, names) = self.parse_header(header)
        self.total_fields = total
        sys.stderr.write("Read "+str(len(names))+\
                             " sample names from VCF header: ")
        sys.stderr.write(str(names)+"\n")

        # next line will be start of VCF body

        # TODO warn when input exceeds eg. 1 TB?
    
        # read VCF body in chunks
        line_count = 0
        chunk_count = 0
        while True:
            lines = infile.readlines(self.buffer_size)
            if lines == []: break
            chunk_count += 1
            for line in lines:
                line_count += 1
                fields = self.parse_body_line(line)
                #stats = update_stats(fields)
        sys.stderr.write("Read "+str(line_count)+" lines from VCF body")
        sys.stderr.write(" in "+str(chunk_count)+" chunk(s).\n")
        return stats

    def parse_body_line(self, line):
        fields = re.split("\s+", line.strip())
        if len(fields) != self.total_fields:
            msg = "Unexpected number of fields in VCF body line; expected "+\
                str(self.total_fields)+", found "+str(len(fields))+\
                " in: "+str(line)
            raise ValueError(msg)
        ref = fields[3]
        alt = fields[4]
        # TODO check the FORMAT field and input to parse_genotype?
        genotypes = [None]*self.total_fields
        #for i in range(self.total_fields):
        #    genotypes[i] = self.parse_genotype(fields[9+i])
        return (ref, alt, genotypes)


    def parse_header(self, column_heads_line):
        # return total headers, and an array of sample names
        fields = re.split("\s+", column_heads_line.strip())
        if len(fields) < 10:
            raise ValueError("No sample names found in column headers: "+\
                                 column_heads_line)
        return (len(fields), fields[9:])

    def parse_genotype(self, input_string):
        # assume format GT:*
        gt_string = re.split(':', input_string).pop(0)
        # TODO get index for pop() from format string
        fields = re.split("\D+", gt_string)
            

# end of vcf_parser

if __name__ == "__main__":
    main()
