Ethnicity calculation using Bayes' Rule
=======================================

Optionally, `vcf_stats.py` estimates the ethnicity of each sample as follows:

- Reads allele frequencies by ethnicity from the INFO field of the VCF body
- Keeps a running total for the log-likelihood of each ethnicity as it reads the VCF file
- Computes the probability of each ethnicity using Bayes' Rule.
