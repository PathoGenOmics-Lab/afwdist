# afwdist

[![PGO badge](https://img.shields.io/badge/PathoGenOmics-Lab-yellow.svg)](https://pathogenomics.github.io/)
[![Release](https://img.shields.io/github/release/PathoGenOmics-Lab/afwdist.svg)](https://github.com/PathoGenOmics-Lab/afwdist/releases)
![Build](https://github.com/PathoGenOmics-Lab/afwdist/actions/workflows/build.yml/badge.svg)
![Test](https://github.com/PathoGenOmics-Lab/afwdist/actions/workflows/test.yml/badge.svg)

An implementation of the pairwise distance metric between groups of genetic variants, based on differences in fixed and non-fixed allele frequencies, described in [Álvarez-Herrera & Sevilla et al. (2024)](https://doi.org/10.1093/ve/veae018) (see also [CITATION.cff](/CITATION.cff)).

Briefly, we define the difference between two vectors of $J$ allele frequencies such that the distance between two samples $M$ and $N$ is the sum for all $I$ polymorphic sites of the differences between the frequency of an allele $j$ at each site $i$:

$$d (M,N) = \sum_{i = 1}^{I} \frac{\sum_{j = 1}^{J} {({{M_{ij}} - {N_{ij}}})}^2} {4 - \sum_{j = 1}^{J} {({{M_{ij}} + {N_{ij}}})}^2}$$

## Usage

### Quick reference

```txt
Usage: afwdist [OPTIONS] --input <INPUT> --reference <REFERENCE> --output <OUTPUT>

Options:
  -i, --input <INPUT>          Input tree in CSV format (mandatory CSV columns are 'sample', 'position', 'sequence' and 'frequency')
  -r, --reference <REFERENCE>  Reference sequence in FASTA format
  -o, --output <OUTPUT>        Output CSV file with distances between each pair of samples
  -v, --verbose                Enable debug messages
  -h, --help                   Print help
  -V, --version                Print version
```

### Inputs and outputs

The software expects a simple table in CSV format (possibly extracted from a VCF file) containing one genetic variants in each row as an input. There must be four columns:

- `sample` (a string): a unique identifier for the group of variants for the pairwise comparison.
- `position` (an integer): the start site of the variant.
- `sequence` (a string): the sequence of the variant (i.e. the alternate allele).
- `frequency` (a real number from 0 to 1): the relative frequency of the variant in the `sample`.

It also expects the reference sequence used for variant calling in FASTA format. As a result, a table in CSV format is produced. This table contains three columns:

- `sample_m` and `sample_n` (strings): the identifiers of the two samples subject to the distance calculation, taken from the input `sample` column.
- `distance` (a real number): the calculated pairwise distance between the samples.

## Citation

> Álvarez-Herrera, M. & Sevilla, J., Ruiz-Rodriguez, P., Vergara, A., Vila, J., Cano-Jiménez, P., González-Candelas, F., Comas, I., & Coscollá, M. (2024). VIPERA: Viral Intra-Patient Evolution Reporting and Analysis. Virus Evolution, 10(1), veae018. https://doi.org/10.1093/ve/veae018
