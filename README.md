oicr_vcf
========

Python code to read a VCF file and output sample statistics in JSON format.

Written for the OICR programming exercise.

Author: Iain Bancarz, iainrb _AT_ gmail _DOT_ com.


Requirements
------------

Requires Python 3.x.


Installation
------------

The file `vcf_stats.py` is self-contained and may be copied to any location of the user's choice.


Usage
-----

Run from the Linux command line as follows:

    ./vcf_stats.py input_file.vcf

Alternatively, use a `-` (dash) to read VCF from standard input:

    ./vcf_stats.py - < input_file.vcf

Output is written to the current directory by default. Specify another
directory with `-o` or `--out`:

    ./vcf_stats.py --out ~/vcf_output input_file.vcf

The optional `-e` or `--ethnicity` argument is a JSON path for output of estimated probabilities of sample ethnicity:

    ./vcf_stats.py --ethnicity ethnicity.json input_file.vcf

See the file `doc/ethnicity.md` for details of the ethnicity estimation.

Run with `-h` or `--help` for more usage information:

    ./vcf_stats.py --help


Output
------

Output is one JSON file for each sample.

File names are of the form `${SAMPLE_NAME}.json`.

Each JSON file contains a single hash with the following keys and values:

- `sample`: The sample name.
- `snps`: Count of each SNP in the input. Indexed by [reference_base] -> [alternate base].
- `ti-tv`: The transition-transversion ratio.
- `variant_count`: Total number of variants found, including SNPs, insertions, deletions, and structural variants.
- `indel_count`: Number of insertions and deletions found.
- `sv_count`: Number of structural variants found.

An example output file appears in `data/example.json`.
