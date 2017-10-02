Ethnicity calculation using Bayes' Theorem
==========================================

Overview
--------

Optionally, `vcf_stats.py` estimates the ethnicity of each sample as follows:

- Reads allele frequencies by ethnicity from the INFO field of the VCF body
- Keeps a running total for the log-likelihood of each ethnicity as it reads the VCF file
- Computes the probability of each ethnicity using Bayes' Theorem.


Bayesian calculation
--------------------

Bayes' Theorem can be used to find the conditional probability of a hypothesis H, given data D:

    Pr(H|D) = (Pr(D|H)Pr(H)) / Pr(D)

where:

- `Pr(H|D)` is the _posterior_ probability of the hypothesis given the data
- `Pr(D|H)` is probability of the data, assuming the hypothesis to be true
- `Pr(H)` is the _prior_ probability of the hypothesis, without the data
- `Pr(D)` is the probability of the data, and acts as a normalizing constant

The INFO column of the VCF file gives allele frequencies by ethnicity, identified by keys `AFR_AF`, `AMR_AF`, `ASN_AF`, `EUR_AF`. This is exactly the probability of observing a variant for a given ethnicity, in other words `Pr(D|H)`. If a given variant is _not_ present, this is also useful data, which occurs with probability `(1 - allele_frequency)`. 

We make the simplifying assumption that all variants are independent, in order to make the computation tractable. In fact, variants are likely to be strongly correlated within an ethnic group; so this is a conservative assumption, likely to produce a conclusion with less certainty than is actually the case.

Given our assumption, we can write:

    Pr(D|H) = Pr(d_1|H) Pr(d_2|H) ... Pr(d_N|H)

where the data points `d_i` are the observations of N variant calls, and Pr(d_i|H) is known from the allele frequency.

Taking logs gives us:

    \log Pr(D|H) = \sum_i \log Pr(d_i|H)

The program `vcf_stats.py` keeps a running total of the log-likelihood for each sample and ethnicity, as it reads the VCF file.

We assume a uniform prior: That is, in the absence of variant call data, all ethnicities are equally likely. So AFR, AMR, ASN and EUR all receive the same prior probability:

    Pr(H) = 0.25

The probability of the data, `Pr(D)` can be computed from the requirement that probabilities `Pr(H|D)` must sum to 1.

We now have enough information to compute the posterior probability of each ethnicity, given the variant calls. The results are written in JSON format.


Note on precision of allele frequency
-------------------------------------

Some alleles are recorded as having frequency 1.0 in a given ethnic group. If interpreted literally, this would mean any sample having one of these alleles is _absolutely certain_ to belong to the ethnic group in question.

However, the literal interpretation is not very realistic. In particular, allele frequencies are recorded with a precision of only two decimal places. So for example, if the 'true' frequency is 0.995, this would be represented as 1.0.

For this reason, allele frequencies of exactly 0 or exactly 1 are represented as 0.995 and 0.005 respectively.


