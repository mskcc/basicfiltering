# basicfiltering

Basic false-positive filters with the following defaults for calls reported by MuTect and VarDict:

1. (Hard Filter) Tumor Total Depth < 5
2. (Hard Filter) Tumor Variant Reads < 3
3. (Hard Filter, VarDict only) Variant call quality < 20 in either tumor or normal
4. (Soft Filter, `f0.01`) Tumor Variant Allele Fraction (VAF) < 1%
5. (Soft Filter, `tnr5`) For non-hotspot loci, ratio of Tumor:Normal VAFs < 5
6. (Soft Filter, `nm2`, VarDict only) Non-hotspot non-SNV with a mean number of mismatches > 2 in tumor
7. (Soft Filter, `mq55`, VarDict only) Non-hotspot non-SNV with a mean mapping quality < 55 in tumor

[![Build Status](https://travis-ci.com/mskcc/basicfiltering.svg?branch=master)](https://travis-ci.com/mskcc/basicfiltering)

### Quick Start 

```bash
# Download and unpack the latest release from GitHub
export GITHUB_URL=`curl -sL https://api.github.com/repos/mskcc/basicfiltering/releases | grep -m1 tarball_url | cut -d\" -f4`
curl -L -o mskcc-basicfiltering.tar.gz $GITHUB_URL; tar -zxf mskcc-basicfiltering.tar.gz; cd mskcc-basicfiltering-*

# Set an environment variable that the cmo package needs
export CMO_RESOURCE_CONFIG=`pwd`/data/cmo_resources.json

# Install dependencies including the mskcc/cmo package needed here
pip install -r requirements.txt

# Read the documentation for required/optional arguments and their default values
python filter_mutect.py --help
python filter_vardict.py --help
python filter_complex.py --help
```
