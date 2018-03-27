# basicfiltering

Basic false-positive filters for VCFs/TXTs from somatic variant callers:

1. Tumor Total Depth >= 5 (default)
2. Tumor Variant Reads >= 3 (default)
3. Tumor Variant Allele Fraction >= 1% (default)
4. Tumor-Normal Variant Allele Fraction (VAF) Ratio >= 5 (default)
5. If known somatic hotspots are provided, the VAF ratio filter above will be skipped for hotspots

[![Build Status](https://travis-ci.org/mskcc/basicfiltering.svg?branch=master)](https://travis-ci.org/mskcc/basicfiltering)

## Requirements:
- pyvcf : [v0.6.8](http://pyvcf.readthedocs.io/en/latest/INTRO.html)
- pandas : [v0.19.2](http://pandas.pydata.org/)
- nose : [v1.3.7](http://nose.readthedocs.io/en/latest/)

## Auto CWL post-process requirements
- Convert inputVcf to have both string and file as input type
- Convert inputTxt to have both string and file as input type
- Convert hotspotVcf to have both string and file as input type

### Works with following versions output formats:

#### MuTect (filter\_mutect.py)
- [MuTect version](https://github.com/broadinstitute/mutect/tree/1.1.4) = 1.1.4
- Takes in the text and vcf file input and filters based on text input

#### PINDEL (filter\_pindel.py)
- [PINDEL version](https://github.com/genome/pindel/tree/v0.2.5a7) = 0.2.5a7
- Takes in a vcf and filters based on it

#### SomaticIndelDetector (filter\_sid.py)
- [SomaticIndelDetector in GATK  version](https://software.broadinstitute.org/gatk/download/) = 2.3-9
- Takes in the text and vcf file input and filters based on text input

#### VarDict (filter\_vardict.py)
- [VarDict version](https://github.com/AstraZeneca-NGS/VarDictJava/tree/v1.4.6) = 1.4.6
- Takes in a vcf and filters based on it
