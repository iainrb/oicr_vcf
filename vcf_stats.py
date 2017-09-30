#! /usr/bin/env python3

# script to parse VCF files, extract basic stats and write in JSON format

# separate JSON file for each individual

import re, sys, json


def main():

    parser = vcf_parser()
    parser.parse_stats(sys.stdin)
    for sample_stats in parser.stats:
        print(json.dumps(sample_stats, sort_keys=True, indent=4))

class vcf_parser:


    def __init__(self):
        self.verbose = True
        self.buffer_size = 1 * 10**6 # input buffer size, in bytes
        self.total_fields = None
        self.total_samples = None
        self.sample_names = []
        self.stats = []

    def find_gt_index(self, format_string):
        # find location of genotype, represented by GT in the format column
        terms = re.split(':', format_string)
        index = None
        i = 0
        for term in terms:
            if term=='GT':
                index = i
                break
            i += 1
        if index == None:
            raise ValueError("Cannot find location of GT in format: "+\
                             format_string)
        return index

    def init_sample_stats(self, name):
        stats = {
            "snps": {
                "A":  {
                    "A": 0,
                    "C": 0,
                    "T": 0,
                    "G": 0,
                    "N": 0
                },
                "C":  {
                    "A": 0,
                    "C": 0,
                    "T": 0,
                    "G": 0,
                    "N": 0
                },
                "T":  {
                    "A": 0,
                    "C": 0,
                    "T": 0,
                    "G": 0,
                    "N": 0
                },

                "G":  {
                    "A": 0,
                    "C": 0,
                    "T": 0,
                    "G": 0,
                    "N": 0
                },
                "N":  {
                    "A": 0,
                    "C": 0,
                    "T": 0,
                    "G": 0,
                    "N": 0
                },
            },
            "sample" : name,
            "ti-tv": 0.0,
            "variant_count": 0,
            "indel_count": 0,
            "sv_count": 0

        }
        return stats
    
    def parse_stats(self, infile):

        # want to read header (as a whole) and body (in chunks)
        # need stats for each individual

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

        (self.total_fields, self.sample_names) = self.parse_header(header)
        self.total_samples = len(self.sample_names)
        
        for i in range(self.total_samples):
            sample_stats = self.init_sample_stats(self.sample_names[i])
            self.stats.append(sample_stats)
            
        if self.verbose:
            sys.stderr.write("Read "+str(len(meta_lines))+\
                             " lines in VCF metadata\n")
            sys.stderr.write("Read "+str(self.total_samples)+\
                             " sample names from VCF header: ")
            sys.stderr.write(str(self.sample_names)+"\n")

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
                (ref, alt, genotypes) = self.parse_body_line(line)
                self.update_stats(ref, alt, genotypes)
        if self.verbose:
            sys.stderr.write("Read "+str(line_count)+" lines from VCF body")
            sys.stderr.write(" in "+str(chunk_count)+" chunk(s).\n")
        return True

    def parse_body_line(self, line):
        fields = re.split("\s+", line.strip())
        if len(fields) != self.total_fields:
            msg = "Unexpected number of fields in VCF body line; expected "+\
                str(self.total_fields)+", found "+str(len(fields))+\
                " in: "+str(line)
            raise ValueError(msg)
        ref = fields[3]
        alt = fields[4]
        gt_index = self.find_gt_index(fields[8])
        # TODO check the FORMAT field and input to parse_genotype?
        genotypes = [None]*self.total_samples
        for i in range(self.total_samples):
            genotypes[i] = self.parse_genotype(fields[9+i], gt_index)
        return (ref, alt, genotypes)

    def parse_header(self, column_heads_line):
        # return total headers, and an array of sample names
        fields = re.split("\s+", column_heads_line.strip())
        if len(fields) < 10:
            raise ValueError("No sample names found in column headers: "+\
                                 column_heads_line)
        return (len(fields), fields[9:])
            
    def parse_genotype(self, input_string, gt_index):
        # legal values: 0,1,. separated by | or /
        permitted_gt = ('0', '1', '.')
        gt_string = re.split(':', input_string)[gt_index]
        genotypes = re.split("[|/]", gt_string)
        if len(genotypes) >= 3:
            raise ValueError("Triploid and higher genotypes not supported")
        for gt in genotypes:
            if gt not in permitted_gt:
                raise ValueError("Illegal genotype character '"+gt+\
                                 "', not in: "+str(permitted_gt))
        return genotypes

    def update_stats(self, ref, alt, genotypes):

        # classify variant as SNP, insertion, deletion, or structural variant
        ref_len = len(ref)
        alt_len = len(alt)
        vartype = None
        if ref_len == 1 and alt_len == 1:
            vartype = 0 # SNP
        elif alt_len > ref_len and ref_len == 1:
            vartype = 1 # insertion
        elif alt_len < ref_len and alt_len == 1:
            vartype = 2 # deletion
        else:
            vartype = 3 # structural variant
        i = 0
        for gt in genotypes:
            for allele_value in gt:
                # ignore '0' for reference, or '.' for no call
                # TODO should we count no calls?
                if allele_value != '1':
                    continue
                elif vartype == 0: # SNP
                    self.stats[i]["snps"][ref][alt] += 1
                    self.stats[i]["variant_count"] += 1
            i += 1

            

# end of class vcf_parser

if __name__ == "__main__":
    main()
