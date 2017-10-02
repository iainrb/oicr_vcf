#! /usr/bin/env python3

# script to parse VCF files, extract basic stats and write in JSON format
# optionally, estimate probability of ethnic groups for each sample

import argparse, math, os, re, sys, json
from decimal import Decimal

def main():
    """Main method to run the VCF stats program"""
    args = construct_argument_parser().parse_args()

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

    if args.ethnicity: enable_ethnicity = True
    else: enable_ethnicity = False

    vcf = vcf_stats(infile, args.verbose, enable_ethnicity)

    if args.infile != '-':
        infile.close()
    
    for sample_stats in vcf.get_sample_stats():
        sample = vcf.get_sample_name(sample_stats)
        outpath = os.path.join(args.out, sample+'.json')
        out = open(outpath, 'w')
        out.write(json.dumps(sample_stats, sort_keys=True, indent=4))
        out.close()
    if args.verbose:
        sys.stderr.write("Wrote JSON output to: "+args.out+"\n")
    if args.ethnicity:
        eth_data = vcf.estimate_ethnicity()
        out = open(args.ethnicity, 'w')
        out.write(json.dumps(eth_data, sort_keys=True, indent=4))
        out.close()
        if args.verbose:
            msg = "Wrote ethnicity output to: "+args.ethnicity+"\n"
            sys.stderr.write(msg)

def construct_argument_parser():
    """ Construct an ArgumentParser object with appropriate options"""
    desc = "Find variant statistics for each sample in a "+\
           "VCF file, and output in JSON format"
    ap = argparse.ArgumentParser(description=desc)
    ap.add_argument('infile',
                    help='Path to input VCF file, or - to read from STDIN')
    ap.add_argument('-o', '--out', metavar='DIR', default=os.getcwd(),
                    help='Directory path for JSON output; defaults to '+\
                    'current working directory. Output filenames are of '+\
                    'the form: ${SAMPLE_NAME}.json')
    ap.add_argument('-e', '--ethnicity', help='Path for output JSON file '+\
                    'containing estimated likelihood of '+\
                    'ethnicities. Optional.')
    ap.add_argument('-v', '--verbose', action='store_true',
                    help='Print additional information to STDERR')
    return ap
    
class vcf_stats:

    """Class to read a VCF file and store statistics and metadata"""

    SNPS_KEY = 'snps'
    VARIANT_COUNT_KEY = 'variant_count'
    INDEL_COUNT_KEY = 'indel_count'
    SV_COUNT_KEY = 'sv_count'
    TI_TV_KEY = 'ti-tv'
    SAMPLE_KEY = 'sample'
    # ethnicity keys
    ASN_KEY = 'ASN'
    AMR_KEY = 'AMR'
    AFR_KEY = 'AFR'
    EUR_KEY = 'EUR'
    
    def __init__(self, infile, verbose, enable_ethnicity):
        """Constructor.

        infile must be a file object; verbose & enable_ethnicity are Boolean.
        Computing log-likelihoods to evaluate the ethnicity requires 
        significant additional runtime; it can be disabled by setting
        enable_ethnicity to False.
        """
        self.verbose = verbose
        self.enable_ethnicity = enable_ethnicity
        self.buffer_size = 50 * 10**6 # input buffer size, in bytes
        self.total_fields = None
        self.total_samples = None
        self.sample_names = []
        self.stats = []
        self.ethnicity_loglik = []
        self.parse_stats(infile)

    def estimate_ethnicity(self):
        """estimate ethnicity of each sample

        uses allele frequencies in the INFO field of the VCF file"""

        # self.ethnicity_loglik = total of log-likelihood by ethnicity
        # use the decimal module for high-precision arithmetic, as with
        # a large number of variants, the log-likelihood is a negative
        # number of large magnitude

        # Apply Bayes' rule to find probability of each ethnicity:
        # H = hypothesis, eg. 'ethnicity is African'
        # D = data, genotype calls and allele frequencies
        # We want probability of hypothesis given data, Pr(H|D)
        #
        # Bayes' rule:  Pr(H|D) = (Pr(D|H)Pr(H)) / Pr(D)
        # where:
        # * Pr(D|H) = exp(self.ethnicity_loglik(H))
        # * Pr(H) = prior probability of ethnicity H (uniform 0.25)
        # * Pr(D) = normalizing constant
        # * Pr(H|D) = (Pr(D|H)Pr(H)) / Pr(D)

        if self.enable_ethnicity == False:
            raise RunTimeError("Ethnicity calculation is not enabled\n")

        prior = {
            self.ASN_KEY: 0.25,
            self.AMR_KEY: 0.25,
            self.AFR_KEY: 0.25,
            self.EUR_KEY: 0.25
        }
        output = []
        for i in range(self.total_samples):
            sample_output = {
                self.ASN_KEY: None,
                self.AMR_KEY: None,
                self.AFR_KEY: None,
                self.EUR_KEY: None
            }
            eth_loglik = self.ethnicity_loglik[i]
            e = Decimal(math.exp(1))
            for key in eth_loglik.keys():
                if key == self.SAMPLE_KEY: continue
                pr_d_h = e**Decimal(eth_loglik[key])
                pr_h = Decimal(prior[key])
                pr_h_d = pr_d_h * pr_h
                msg = "\t"+str(pr_d_h)+"\n"
                sample_output[key] = pr_h_d
            # normalize probabilities so they sum to 1
            total = sum(sample_output.values())
            for key in sample_output.keys():
                try:
                    normalized = sample_output[key] / total
                    sample_output[key] = float(normalized)
                except ZeroDivisionError:
                    # eg. INFO field had no ethnic allele frequencies
                    sample_output[key] = 0.0
            sample_output[self.SAMPLE_KEY] = self.sample_names[i]
            output.append(sample_output)
        return output
        
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

    def get_sample_name(self, sample_stats):
        """Get the sample name, from that sample's stats list entry"""
        return sample_stats[self.SAMPLE_KEY]

    def get_sample_stats(self):
        """Return the stats list, for subsequent processing and output"""
        return self.stats

    def init_ethnicity_loglik(self, name):
        """initialise an empty data structure with the sample name

        dictionary contains log-likelihood running totals by ethnicity
        """
        stats = {
            self.SAMPLE_KEY: name,
            self.AMR_KEY: 0.0,
            self.ASN_KEY: 0.0,
            self.AFR_KEY: 0.0,
            self.EUR_KEY: 0.0
        }
        return stats
    
    def init_sample_stats(self, name):
        """initialise an empty data structure with the sample name"""
        stats = {
            self.SNPS_KEY: {
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
            self.SAMPLE_KEY: name,
            self.TI_TV_KEY: 0.0,
            self.VARIANT_COUNT_KEY: 0,
            self.INDEL_COUNT_KEY: 0,
            self.SV_COUNT_KEY: 0
        }
        return stats

    def is_transition(self, ref, alt):
        """Is the given SNP an transition?

        Return status of the given reference and alternate alleles:
        - transition: A->G, G->A, C->T, T->C
        - transversion: A->C, C->A, A->T, T->A, G->T, T->G, G->C, C->G"""
        status = False
        if (ref=='A' and alt=='G') or (ref=='G' and alt=='A') or \
           (ref=='C' and alt=='T') or (ref=='T' and alt=='C'):
            status = True
        return status
    
    def parse_stats(self, infile):
        """Read a VCF file and populate instance variables"""
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
            if self.enable_ethnicity:
                eth_stats = self.init_ethnicity_loglik(self.sample_names[i])
                self.ethnicity_loglik.append(eth_stats)
        if self.verbose:
            sys.stderr.write("Read "+str(len(meta_lines))+\
                             " lines in VCF metadata\n")
            sys.stderr.write("Read "+str(self.total_samples)+\
                             " sample names from VCF header\n")
        # read VCF body in chunks
        line_count = 0
        block_count = 0
        while True:
            lines = infile.readlines(self.buffer_size)
            if lines == []: break
            block_count += 1
            for line in lines:
                line_count += 1
                (ref, alts, genotypes, info) = self.parse_body_line(line)
                self.update_stats(ref, alts, genotypes, info)
        if self.verbose:
            sys.stderr.write("Read "+str(line_count)+" lines from VCF body"+\
                             " in "+str(block_count)+" block(s) of maximum "+\
                             str(self.buffer_size)+" bytes\n")
        # update with transition/transversion ratios
        self.update_sample_titv()     
        return True

    def parse_body_line(self, line):
        """parse a line from the body of a VCF file

        returns: reference, one or more alternates, genotypes, info
        info is a dictionary of ethnic allele frequences from the INFO field
        """
        fields = re.split("\s+", line.strip())
        if len(fields) != self.total_fields:
            msg = "Unexpected number of fields in VCF body line; expected "+\
                str(self.total_fields)+", found "+str(len(fields))+\
                " in: "+str(line)
            raise VCFInputError(msg)
        ref = fields[3]
        alts = re.split(',', fields[4])
        if self.enable_ethnicity:
            info = self.parse_info(fields[7])
        else:
            info = {}
        gt_index = self.find_gt_index(fields[8])
        genotypes = [None]*self.total_samples
        for i in range(self.total_samples):
            genotypes[i] = self.parse_genotype(fields[9+i], gt_index)
        return (ref, alts, genotypes, info)

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

        Sample field consists of one or more colon-delimited sub-fields.
        Legal values for the genotype sub-field:
        0,1,. separated by | or /"""
        permitted_gt = ('0', '1', '.')
        gt_string = re.split(':', input_string)[gt_index]
        genotypes = re.split("[|/]", gt_string)
        if len(genotypes) >= 3:
            raise VCFInputError("Polyploid genotypes not supported")
        for gt in genotypes:
            if re.match('<.*>', gt):
                raise VCFInputError("ID string for alternate not supported")
            elif re.search('\[|\]', gt):
                raise VCFInputError("Breakends for alternate not supported")
            elif re.search('[^0-9\.]', gt):
                raise VCFInputError("Illegal genotype character in '"+gt+\
                                    "', not an integer or '.'")
        return genotypes

    def parse_info(self, info_string):
        """Parse the INFO field in VCF body; get allele frequencies"""
        fields = re.split(';', info_string)
        af_keys = ('AMR_AF', 'ASN_AF', 'AFR_AF', 'EUR_AF')
        info = {}
        for field in fields:
            try:
                (key, value) = re.split('=', field)
                if key in af_keys:
                    info[key] = float(value)
            except ValueError: # eg. no = sign in field
                continue
        return info
    
    def update_sample_titv(self):
        """Update transition-transversion ratio for each sample"""
        for i in range(self.total_samples):
            titv_ratio = self.titv(self.stats[i])
            self.stats[i][self.TI_TV_KEY] = titv_ratio

    def update_ethnicity(self, sample_index, is_variant, info):
        """Update running totals for log-likelihood of ethnicity"""
        key_map = {
            'AMR_AF': self.AMR_KEY,
            'ASN_AF': self.ASN_KEY,
            'AFR_AF': self.AFR_KEY,
            'EUR_AF': self.EUR_KEY,
        }
        for key in info.keys():
            # correction for very high/low allele frequency
            # VCF INFO field is only precise to 2 d.p.
            # so 1.0 has the same representation as 0.995
            af = info[key]
            delta = 0.0001
            if math.fabs(1.0 - af) < delta: af = 0.995
            elif math.fabs(af) < delta: af = 0.005
            # update running total for variant or no-variant
            if is_variant:
                loglik = math.log(af) # ln Pr(variant|ethnicity)
            else:
                loglik = math.log(1.0 - af) # ln Pr(no-variant|ethnicity)
            self.ethnicity_loglik[sample_index][key_map[key]] += loglik
            
    def update_stats(self, ref, alts, genotypes, info):
        """Update running totals for all samples, for a given VCF line

        Also updates ethnicity log-likelihood totals, using 'info' argument
        """
        # classify alterantes as SNP, indel, or structural variant
        ref_len = len(ref)
        alt_types = []
        for alt in alts:
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
            alt_types.append(vartype)
        i = 0
        for gt in genotypes:
            # find which type of variant is present
            # does not support multiple variant types in same genotype,
            # eg. SNP on one chromosome and indel on the other
            variants = set()
            for allele_value in gt: # for each chromosome
                if self.update_ethnicity:
                    if allele_value == '.':
                        pass
                    elif allele_value == '0':
                        self.update_ethnicity(i, False, info)
                    else:
                        self.update_ethnicity(i, True, info)
                # ignore '0' for reference, or '.' for no call
                if allele_value != '.' and allele_value != '0':
                    alt_index = int(allele_value) - 1
                    variants.add(alt_index)
                    if alt_types[alt_index] == 0: # SNP
                        alt = alts[alt_index]
                        self.stats[i][self.SNPS_KEY][ref][alt] += 1
            # each variant type is counted only once in the variant total
            # eg. homozygous SNP (two alternate alleles) is not double-counted
            # however a SNP and indel at same position are counted separately
            for variant_index in variants:
                vartype = alt_types[variant_index]
                self.stats[i][self.VARIANT_COUNT_KEY] += 1
                if vartype == 1 or vartype == 2:
                     self.stats[i][self.INDEL_COUNT_KEY] += 1
                elif vartype == 3:
                     self.stats[i][self.SV_COUNT_KEY] += 1
            i += 1

    def titv(self, sample_stats):
        """find the ti-tv (transition-transversion) ratio"""
        ti = 0
        tv = 0
        counts = sample_stats[self.SNPS_KEY]
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

# end of class vcf_parser

    
class VCFInputError(Exception):
    """Error for badly formed VCF input"""

    pass


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
