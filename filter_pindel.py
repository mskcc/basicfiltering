#!/usr/bin/env python
'''
@description : This tool helps to filter a pindel vcf
@created : 07/17/2016
@author : Ronak H Shah, Cyriac Kandoth, Zuojian Tang

'''

from __future__ import division
import argparse, sys, os, time, logging, cmo

logging.basicConfig(
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logging.DEBUG)
logger = logging.getLogger('filter_pindel')
try:
    import vcf
    from vcf.parser import _Info as VcfInfo, _Format as VcfFormat
except ImportError:
    logger.fatal("filter_mutect: pyvcf is not installed")
    sys.exit(1)

def main():
    parser = argparse.ArgumentParser(prog='filter_pindel.py',description='Filter indels from the output of pindel',usage='%(prog)s [options]')
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",help="make lots of noise")
    parser.add_argument("-ivcf", "--inputVcf", action="store", dest="inputVcf", required=True, type=str, metavar='SomeID.vcf',help="Input vcf file which needs to be filtered")
    parser.add_argument("-tsn", "--tsampleName", action="store", dest="tsampleName", required=True, type=str, metavar='SomeName',help="Name of the tumor Sample")
    parser.add_argument("-rf", "--refFasta", action="store", dest="refFasta", required=True, type=str, metavar='ref.fa', help="Reference genome in fasta format")
    parser.add_argument("-dp", "--totaldepth", action="store", dest="dp", required=False, type=int, default=5, metavar='5',help="Tumor total depth threshold")
    parser.add_argument("-ad", "--alleledepth", action="store", dest="ad", required=False, type=int, default=3, metavar='3',help="Tumor allele depth threshold")
    parser.add_argument("-tnr", "--tnRatio", action="store", dest="tnr", required=False, type=int, default=5, metavar='5',help="Tumor-Normal variant fraction ratio threshold ")
    parser.add_argument("-vf", "--variantfraction", action="store", dest="vf", required=False, type=float, default=0.01, metavar='0.01',help="Tumor variant fraction threshold ")
    parser.add_argument("-min", "--min_var_len", action="store", dest="min", required=False, default=0, type=int, metavar='0',help="Minimum length of the indels")
    parser.add_argument("-max", "--max_var_len", action="store", dest="max", required=False, default=200, type=int, metavar='200',help="Max length of the indels")
    parser.add_argument("-mhl", "--max_hom_len", action="store", dest="mhl", required=False, default=5, type=int, metavar='5',help="Max length of micro-homology at indel breakpoint")
    parser.add_argument("-hvcf", "--hotspotVcf", action="store", dest="hotspotVcf", required=False, type=str, metavar='hotspot.vcf',help="Input vcf file with hotspots that skip VAF ratio filter")
    parser.add_argument("-o", "--outDir", action="store", dest="outdir", required=False, type=str, metavar='/somepath/output',help="Full Path to the output dir.")

    args = parser.parse_args()
    if(args.verbose):
        logger.info("Started the run for doing standard filter.")
    RunStdFilter(args)
    if(args.verbose):
        logger.info("Finished the run for doing standard filter.")

def RunStdFilter(args):
    vcf_out = os.path.basename(args.inputVcf)
    vcf_out = os.path.splitext(vcf_out)[0]
    if(args.outdir):
        vcf_out = os.path.join(args.outdir,vcf_out)
    txt_out = vcf_out + "_STDfilter.txt"
    vcf_out = vcf_out + "_STDfilter.vcf"
    vcf_reader = vcf.Reader(open(args.inputVcf, 'r'))
    del vcf_reader.infos['END']
    vcf_reader.infos['set'] = VcfInfo('set', '.', 'String', 'The variant callers that reported this event', 'mskcc/basicfiltering', 'v0.2.2')
    vcf_reader.formats['DP'] = VcfFormat('DP', '1', 'Integer', 'Total read depth at this site')
    vcf_reader.formats['AD'] = VcfFormat('AD', 'R', 'Integer', 'Allelic depths for the ref and alt alleles in the order listed')
    vcf_reader.formats['PL'] = VcfFormat('PL', 'G', 'Integer', 'Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification')

    allsamples = list(vcf_reader.samples)
    if(len(allsamples) != 2):
        if(args.verbose):
            logger.critical("The VCF does not have two genotype columns. Please input a proper vcf with Tumor/Normal columns")
        sys.exit(1)

    # If the caller reported the normal genotype column before the tumor, swap those around
    if_swap_sample = False
    if(allsamples[1] == args.tsampleName):
        if_swap_sample = True
        vcf_reader.samples[0] = allsamples[1]
        vcf_reader.samples[1] = allsamples[0]
    nsampleName = vcf_reader.samples[1]

    vcf_writer = vcf.Writer(open(vcf_out, 'w'), vcf_reader)
    txt_fh = open(txt_out, "wb")

    # If provided, load hotspots into a dictionary for quick lookup
    hotspot = {}
    if(args.hotspotVcf):
        hvcf_reader = vcf.Reader(open(args.hotspotVcf, 'r'))
        for record in hvcf_reader:
            genomic_locus = str(record.CHROM) + ":" + str(record.POS)
            hotspot[genomic_locus] = True

    for record in vcf_reader:
        tcall = record.genotype(args.tsampleName)
        recordType = record.INFO['SVTYPE']
        recordLen = abs(int(record.INFO['SVLEN']))
        recordHomLen = int(record.INFO['HOMLEN'])
        record.INFO.pop('END',0)
        locus = str(record.CHROM) + ":" + str(record.POS)
        # Remove INFO/SVTYPE altogether because it causes VEP to avoid properly annotating it
        del record.INFO['SVTYPE']

        if(tcall['AD'] is not None):
            trd, tad = tcall['AD']
        else:
            trd = 0
            tad = 0
        tdp = trd + tad
        tvf = int(tad)/float(tdp) if(tdp != 0) else 0
        # print "Tdata: ",trd,tad
        ncall = record.genotype(nsampleName)
        if(ncall):
            if(ncall['AD'] is not None):
                nrd, nad = ncall['AD']
            else:
                nrd = 0
                nad = 0
            ndp = nrd + nad
            nvf = nad/ndp if(ndp != 0) else 0
            nvfRF = int(args.tnr) * nvf

        if (if_swap_sample):
            nrm = record.samples[0]
            tum = record.samples[1]
            record.samples[0] = tum
            record.samples[1] = nrm
        if((recordLen >= int(args.min)) and (recordLen <= int(args.max)) and (recordHomLen <= int(args.mhl))):
            record.add_info('set', 'Pindel')
            if tvf > nvfRF or locus in hotspot:
                if((tdp >= int(args.dp)) & (tad >= int(args.ad)) & (tvf >= float(args.vf))):
                    vcf_writer.write_record(record)
                    vcf_writer.flush()
                    txt_fh.write( args.tsampleName + "\t" + record.CHROM + "\t" + str(record.POS) + "\t" + str(record.REF) + "\t" + str(record.ALT[0]) + "\t" + "." + "\n")
    vcf_writer.close()
    txt_fh.close()
    # Normalize the events in the VCF, produce a bgzipped VCF, then tabix index it
    norm_gz_vcf = cmo.util.normalize_vcf(vcf_out, args.refFasta)
    cmo.util.tabix_file(norm_gz_vcf)
    return(norm_gz_vcf)

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    logging.info("filter_pindel: Elapsed time was %g seconds", totaltime)
    sys.exit(0)
