#! /usr/bin/env python3

# script to parse VCF files, extract basic stats and write in JSON format

import argparse, re, os, sys, json

def main():
    """Main method to run the VCF stats program"""
    desc = "Find variant statistics for each sample in a "+\
           "VCF file, and output in JSON format"
    ap = argparse.ArgumentParser(description=desc)
    ap.add_argument('infile',
                    help='Path to input VCF file, or - to read from STDIN')
    ap.add_argument('-o', '--out', metavar='DIR', default=os.getcwd(),
                    help='Directory path for JSON output; defaults to '+\
                    'current working directory. Output filenames are of '+\
                    'the form: ${SAMPLE_NAME}.json')
    ap.add_argument('-v', '--verbose', action='store_true',
                    help='Print additional information to STDERR')
    args = ap.parse_args()

    infile = None
    if not os.path.exists(args.out):
        raise ValueError("Output path '"+args.out+"' does not exist")
    elif not os.path.isdir(args.out):
        raise ValueError("Output path '"+args.out+"' is not a directory")
    if args.infile == '-':
        infile = sys.stdin
    elif not os.path.exists(args.infile):
        raise ValueError("Input path '"+args.infile+"' does not exist")
    elif not os.path.isfile(args.infile):
        raise ValueError("Input path '"+args.infile+"' is not a regular file")
    else:
        infile = open(args.infile, 'r')

    vcf = vcf_stats(infile, args.verbose)

    if args.infile != '-':
        infile.close()
    
    for sample_stats in vcf.stats:
        sample = sample_stats['sample']
        outpath = os.path.join(args.out, sample+'.json')
        out = open(outpath, 'w')
        out.write(json.dumps(sample_stats, sort_keys=True, indent=4))
        out.close()


class vcf_stats:

    """Class to read a VCF file and contain statistics and metadata"""

    def __init__(self, infile, verbose):
        """Constructor. infile must be a file object; verbose is Boolean"""
        self.verbose = verbose
        self.buffer_size = 1 * 10**6 # input buffer size, in bytes
        self.total_fields = None
        self.total_samples = None
        self.sample_names = []
        self.stats = []
        self.parse_stats(infile)

    def find_gt_index(self, format_string):
        """parse the VCF format field, to find location of the genotype"""
        terms = re.split(':', format_string)
        index = None
        i = 0
        for term in terms:
            if term=='GT':
                index = i
                break
            i += 1
        if index == None:
            raise VCFInputError("Cannot find location of GT in format: "+\
                             format_string)
        return index

    def init_sample_stats(self, name):
        """initialise an empty data structure with the sample name"""
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

    def is_transition(self, ref, alt):
        """Is the given SNP an transition?
        - transition: A->G, G->A, C->T, T->C
        - transversion: A->C, C->A, A->T, T->A, G->T, T->G, G->C, C->G"""
        status = False
        if (ref=='A' and alt=='G') or (ref=='G' and alt=='A') or \
           (ref=='C' and alt=='T') or (ref=='T' and alt=='C'):
            status = True
        return status
    
    def parse_stats(self, infile):
        """Read metadata, header, and body of a VCF file and populate
        instance variables"""
        # read VCF meta lines and header
        meta_lines = []
        header = None
        while True:
            # read header lines one at a time
            line = infile.readline()
            if line == '':
                msg = "Reached end of file without finding end of VCF header"
                raise VCFInputError(msg)
            if re.match('##', line):
                meta_lines.append(line)
            elif re.match('#CHROM', line):
                header = line
                break
            else:
                msg = "Unexpected line in VCF header; line "+\
                    "does not start with '##' or '#CHROM': "+line
                raise VCFInputError(msg)
        (self.total_fields, self.sample_names) = self.parse_header(header)
        self.total_samples = len(self.sample_names)
        for i in range(self.total_samples):
            sample_stats = self.init_sample_stats(self.sample_names[i])
            self.stats.append(sample_stats)
        if self.verbose:
            sys.stderr.write("Read "+str(len(meta_lines))+\
                             " lines in VCF metadata\n")
            sys.stderr.write("Read "+str(self.total_samples)+\
                             " sample names from VCF header\n")
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

        # update with transition/transversion ratios
        self.update_sample_titv()        
        return True

    def parse_body_line(self, line):
        """parse a line from the body of a VCF file"""
        fields = re.split("\s+", line.strip())
        if len(fields) != self.total_fields:
            msg = "Unexpected number of fields in VCF body line; expected "+\
                str(self.total_fields)+", found "+str(len(fields))+\
                " in: "+str(line)
            raise VCFInputError(msg)
        ref = fields[3]
        alt = fields[4]
        gt_index = self.find_gt_index(fields[8])
        # TODO check the FORMAT field and input to parse_genotype?
        genotypes = [None]*self.total_samples
        for i in range(self.total_samples):
            genotypes[i] = self.parse_genotype(fields[9+i], gt_index)
        return (ref, alt, genotypes)

    def parse_header(self, column_heads_line):
        """Parse the header line of a VCF file.
        Return total headers, and an array of sample names."""
        fields = re.split("\s+", column_heads_line.strip())
        if len(fields) < 10:
            raise VCFInputError("No sample names found in column headers: "+\
                                 column_heads_line)
        names = fields[9:]
        name_set = set()
        for name in names:
            if name in name_set:
                raise VCFInputError("Sample name '"+name+\
                                    "' appears more than once in VCF header")
            else:
                name_set.add(name)
        return (len(fields), names)
            
    def parse_genotype(self, input_string, gt_index):
        """Find genotype from a sample field in a VCF file.
        Legal values for the genotype sub-field:
        0,1,. separated by | or /"""
        permitted_gt = ('0', '1', '.')
        gt_string = re.split(':', input_string)[gt_index]
        genotypes = re.split("[|/]", gt_string)
        if len(genotypes) >= 3:
            raise VCFInputError("Polyploid genotypes not supported")
        for gt in genotypes:
            if gt not in permitted_gt:
                raise VCFInputError("Illegal genotype character '"+gt+\
                                    "', not in: "+str(permitted_gt))
        return genotypes

    def update_sample_titv(self):
        """Update transition-transversion ratio for each sample"""
        for i in range(self.total_samples):
            titv_ratio = self.titv(self.stats[i])
            self.stats[i]["ti-tv"] = titv_ratio
    
    def update_stats(self, ref, alt, genotypes):
        """Update running totals for all samples, for a given VCF line"""
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
            is_variant = False
            for allele_value in gt:
                # ignore '0' for reference, or '.' for no call
                if allele_value != '1':
                    continue
                else:
                    is_variant = True
                if vartype == 0: # SNP
                    self.stats[i]["snps"][ref][alt] += 1
            # variant may be homozygous or heterozygous
            # each type is counted only once in variant total
            if is_variant:
                self.stats[i]["variant_count"] += 1
                if vartype == 1 or vartype == 2:
                     self.stats[i]["indel_count"] += 1
                elif vartype == 3:
                     self.stats[i]["sv_count"] += 1
            i += 1

    def titv(self, sample_stats):
        """find the ti-tv (transition-transversion) ratio"""
        ti = 0
        tv = 0
        counts = sample_stats["snps"]
        bases = ('A', 'C', 'G', 'T')
        for ref in bases:
            for alt in bases:
                if self.is_transition(ref, alt): ti += counts[ref][alt]
                else: tv += counts[ref][alt]
        titv = None
        try:
            titv = float(ti) / tv
        except ZeroDivisionError:
            titv = 0.0
        return titv


class VCFInputError(Exception):
    """Error for badly formed VCF input"""
    pass

# end of class vcf_parser

if __name__ == "__main__":
    main()

"""Author: Iain Bancarz <iainrb _AT_ gmail _DOT_ com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
