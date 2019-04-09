# basicfiltering

Basic false-positive filters for VCFs/TXTs from somatic variant callers:

1. Tumor Total Depth < 5
2. Tumor Variant Reads < 3
3. Tumor Variant Allele Fraction (VAF) < 1%
4. For non-hotspot loci, ratio of Tumor:Normal VAFs < 5

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
