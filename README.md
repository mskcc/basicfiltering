# basicfiltering
Basic Filtering for: 

1. 1% Variant Allele Frequency and 
2. 3 Variant Reads and 
3. Tumor-Normal Variant Allele Frequency Ratio to be greater then 5 times 

for Multiple Tools

## Works with following versions output formats:

#### SomaticIndelDetector (filter\_sid.py)
- [SomaticIndelDetector in GATK  version](https://software.broadinstitute.org/gatk/download/) = 2.3-9
- Takes in the text and vcf file input and filters based on text input. 	

```
usage: filter_sid.py [options]

Filter Indels from the output of SomaticIndelDetector

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         make lots of noise [default]
  -ivcf SomeID.vcf, -inputVcf SomeID.vcf
                        Input SomaticIndelDetector vcf file which needs to be
                        filtered
  -itxt SomeID.txt, -inputTxt SomeID.txt
                        Input SomaticIndelDetector txt file which needs to be
                        filtered
  -tsn SomeName, --tsampleName SomeName
                        Name of the tumor Sample
  -dp 0, --totaldepth 0
                        Tumor total depth threshold
  -ad 5, --alleledepth 5
                        Tumor allele depth threshold
  -tnr 5, --tnRatio 5   Tumor-Normal variant frequency ratio threshold
  -vf 0.01, --variantfrequency 0.01
                        Tumor variant frequency threshold
  -hvcf hostpot.vcf, --hotspotVcf hostpot.vcf
                        Input bgzip / tabix indexed hotspot vcf file to used
                        for filtering
  -o /somepath/output, --outDir /somepath/output
                        Full Path to the output dir.

```

#### MuTect (filter\_mutect.py)
- [MuTect version](https://github.com/broadinstitute/mutect/tree/1.1.4) = 1.1.4
- Takes in the text and vcf file input and filters based on text input.

```
usage: filter_mutect.py [options]

Filter SNPS from the output of muTect v1.14

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         make lots of noise [default]
  -ivcf SomeID.vcf, -inputVcf SomeID.vcf
                        Input vcf muTect file which needs to be filtered
  -itxt SomeID.txt, -inputTxt SomeID.txt
                        Input txt muTect file which needs to be filtered
  -tsn SomeName, --tsampleName SomeName
                        Name of the tumor Sample
  -dp 0, --totaldepth 0
                        Tumor total depth threshold
  -ad 5, --alleledepth 5
                        Tumor allele depth threshold
  -tnr 5, --tnRatio 5   Tumor-Normal variant frequency ratio threshold
  -vf 0.01, --variantfrequency 0.01
                        Tumor variant frequency threshold
  -hvcf hostpot.vcf, --hotspotVcf hostpot.vcf
                        Input bgzip / tabix indexed hotspot vcf file to used
                        for filtering
  -o /somepath/output, --outDir /somepath/output
                        Full Path to the output dir.

```

#### VarDict (filter\_vardict.py)
- [VarDict version](https://github.com/AstraZeneca-NGS/VarDictJava/tree/v1.4.6) = 1.4.6
- Takes in a vcf and filters based on it

```
usage: filter_vardict.py [options]

Filter Indels from the output of vardict

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         make lots of noise [default]
  -i SomeID.vcf, -inputVcf SomeID.vcf
                        Input vcf vardict file which needs to be filtered
  -tsn SomeName, --tsampleName SomeName
                        Name of the tumor Sample
  -dp 0, --totaldepth 0
                        Tumor total depth threshold
  -ad 5, --alleledepth 5
                        Tumor allele depth threshold
  -tnr 5, --tnRatio 5   Tumor-Normal variant frequency ratio threshold
  -vf 0.01, --variantfrequency 0.01
                        Tumor variant frequency threshold
  -hvcf hostpot.vcf, --hotspotVcf hostpot.vcf
                        Input bgzip / tabix indexed hotspot vcf file to used
                        for filtering
  -o /somepath/output, --outDir /somepath/output
                        Full Path to the output dir.
```
#### PINDEL (filter\_pindel.py)
- [PINDEL version](https://github.com/genome/pindel/tree/v0.2.5a7) = 0.2.5a7
- Takes in a vcf and filters based on it

```
usage: filter_pindel.py [options]

Filter Indels from the output of pindel

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         make lots of noise [default]
  -i SomeID.vcf, -inputVcf SomeID.vcf
                        Input vcf freebayes file which needs to be filtered
  -tsn SomeName, --tsampleName SomeName
                        Name of the tumor Sample
  -dp 0, --totaldepth 0
                        Tumor total depth threshold
  -ad 5, --alleledepth 5
                        Tumor allele depth threshold
  -tnr 5, --tnRatio 5   Tumor-Normal variant frequency ratio threshold
  -vf 0.01, --variantfrequency 0.01
                        Tumor variant frequency threshold
  -o /somepath/output, --outDir /somepath/output
                        Full Path to the output dir.
  -min 25, --min_var_len 25
                        Minimum length of the Indels
  -max 500, --max_var_len 500
                        Max length of the Indels
  -hvcf hostpot.vcf, --hotspotVcf hostpot.vcf
                        Input bgzip / tabix indexed hotspot vcf file to used
                        for filtering

```

