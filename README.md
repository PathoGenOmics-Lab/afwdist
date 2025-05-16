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

The program takes as input a table in CSV format (possibly derived from a VCF file) where each row represents a single genetic variant. The input table must contain four columns:

- `sample` (a string): a unique identifier for the group of variants used in pairwise comparisons.
- `position` (an integer): the site of the variant.
- `sequence` (a string): the sequence of the variant (i.e. the alternate allele).
- `frequency` (a real number from 0 to 1): the relative frequency of the variant within the sample.

In addition to the variant table, the program requires a reference sequence in FASTA format. The sequence should be the same one used for variant calling. This reference is used to infer the frequencies of reference alleles, assuming that any frequency not taken up by listed variants belongs to the reference allele at that site. In addition to the pairwise distance between samples, the distance between each sample and the reference sequence is also calculated by building a reference sample as a baseline with no variant alleles (i.e. all sites are assumed to have an allele frequency of 1).

The distance of each sample is calculated against the reference as well, treating it as a normal sample with no allele variants (all reference allele frequencies are fixed within the reference virtual sample).

As a result, a table in CSV format is produced. This table contains three columns:

- `sample_m` and `sample_n` (strings): the identifiers of the two samples being compared.
- `distance` (a real number): the calculated pairwise distance between the two samples.

## Citation

> Álvarez-Herrera, M. & Sevilla, J., Ruiz-Rodriguez, P., Vergara, A., Vila, J., Cano-Jiménez, P., González-Candelas, F., Comas, I., & Coscollá, M. (2024). VIPERA: Viral Intra-Patient Evolution Reporting and Analysis. Virus Evolution, 10(1), veae018. https://doi.org/10.1093/ve/veae018

## Contributors

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tbody>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/ahmig"><img src="https://avatars.githubusercontent.com/u/30174538?v=4?s=100" width="100px;" alt="Miguel Álvarez Herrera"/><br /><sub><b>Miguel Álvarez Herrera</b></sub></a><br /><a href="https://github.com/ahmig/afwdist/commits?author=ahmig" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/SeviJordi"><img src="https://avatars.githubusercontent.com/u/124465736?v=4?s=100" width="100px;" alt="Jordi Sevilla Fortuny"/><br /><sub><b>Jordi Sevilla Fortuny</b></sub></a><br /><a href="https://github.com/ahmig/afwdist/issues?q=author%3ASeviJordi" title="Bug reports">🐛</a> <a href="#userTesting-SeviJordi" title="User Testing">📓</a></td>
    </tr>
  </tbody>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification.
