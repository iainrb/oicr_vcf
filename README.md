oicr_vcf
========

Python code to read a VCF file and output sample statistics in JSON format.

Written for the OICR programming exercise.

Author: Iain Bancarz, iainrb _AT_ gmail _DOT_ com.


Requirements
------------

Requires Python 3.x.


Usage
-----

Run from the Linux command line as follows:

    ./vcf_stats.py input_file.vcf

Alternatively, use a `-` (dash) to read VCF from standard input:

    ./vcf_stats.py - < input_file.vcf

Output is written to the current directory by default. Specify another
directory with `-o` or `--out`:

    ./vcf_stats.py --out ~/vcf_output input_file.vcf

Run with `-h` or `--help` for more information:

    ./vcf_stats.py --help


Output
------

Output is one JSON file for each sample.

File names are of the form `${SAMPLE_NAME}.json`.

Each JSON file contains a single hash with the following keys and values:

- `sample`: The sample name.
- `snps`: Count of each SNP in the input. Indexed by [reference_base] -> [alternate base]. Both alleles in a heterozygous SNP will be counted; for example AA->GC counts as one A->G, and one A->C. For a heterozygous SNP, only the variant allele is counted, so AA->AC is recorded as one A->C.
- `ti-tv`: The transition-transversion ratio.
- `variant count`: Total number of variants found, including SNPs, insertions, deletions, and structural variants. (A heterozygous SNP counts as only one variant.)
- `indel_count`: Number of insertions and deletions found.
- `sv_count`: Number of structural variants found.