#!/usr/bin/env python
'''
@Description : This tool helps to filter SomaticIndelDetector (GATK version: 2.3-9) vcf and txt
@Created : 07/25/2016
@Updated : 03/26/2018
@author : Ronak H Shah, Cyriac Kandoth

'''
from __future__ import division
import argparse, sys, os, time, logging

logging.basicConfig(
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logging.DEBUG)
logger = logging.getLogger('filter_sid')
try:
    import vcf
except ImportError:
    logger.fatal("filter_sid: pyvcf is not installed, please install pyvcf as it is required to run the mapping.")
    sys.exit(1)

def main():
    parser = argparse.ArgumentParser(prog='filter_sid.py',description='Filter indels from the output of SomaticIndelDetector in GATK v2.3-9',usage='%(prog)s [options]')
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",help="make lots of noise")
    parser.add_argument("-ivcf", "--inputVcf", action="store", dest="inputVcf", required=True, type=str, metavar='SomeID.vcf',help="Input SomaticIndelDetector vcf file which needs to be filtered")
    parser.add_argument("-itxt", "--inputTxt", action="store", dest="inputTxt", required=True, type=str, metavar='SomeID.txt',help="Input SomaticIndelDetector txt file which needs to be filtered")
    parser.add_argument("-tsn", "--tsampleName", action="store", dest="tsampleName", required=True, type=str, metavar='SomeName',help="Name of the tumor Sample")
    parser.add_argument("-dp", "--totaldepth", action="store", dest="dp", required=False, type=int, default=5, metavar='5',help="Tumor total depth threshold")
    parser.add_argument("-ad", "--alleledepth", action="store", dest="ad", required=False, type=int, default=3, metavar='3',help="Tumor allele depth threshold")
    parser.add_argument("-tnr", "--tnRatio", action="store", dest="tnr", required=False, type=int, default=5, metavar='5',help="Tumor-Normal variant fraction ratio threshold")
    parser.add_argument("-vf", "--variantfraction", action="store", dest="vf", required=False, type=float, default=0.01, metavar='0.01',help="Tumor variant fraction threshold")
    parser.add_argument("-hvcf", "--hotspotVcf", action="store", dest="hotspotVcf", required=False, type=str, metavar='hotspot.vcf',help="Input vcf file with hotspots that skip VAF ratio filter")
    parser.add_argument("-o", "--outDir", action="store", dest="outdir", required=False, type=str, metavar='/somepath/output',help="Full Path to the output dir.")

    args = parser.parse_args()
    if(args.verbose):
        logger.info("Started the run for doing standard filter.")
    (stdfilterVCF) = RunStdFilter(args)
    if(args.verbose):
        logger.info("Finished the run for doing standard filter.")

def RunStdFilter(args):
    vcf_out = os.path.basename(args.inputVcf)
    vcf_out = os.path.splitext(vcf_out)[0]
    txt_out = vcf_out
    if(args.outdir):
        vcf_out = os.path.join(args.outdir,vcf_out + "_STDfilter.vcf")
        txt_out = os.path.join(args.outdir,txt_out + "_STDfilter.txt")
    else:
        vcf_out = vcf_out + "_STDfilter.vcf"
        txt_out = txt_out + "_STDfilter.txt"
    vcf_reader = vcf.Reader(open(args.inputVcf, 'r'))
    vcf_writer = vcf.Writer(open(vcf_out, 'w'), vcf_reader)
    txt_fh = open(txt_out, "wb")
    allsamples = vcf_reader.samples
    if(len(allsamples) == 2):
        sample1 = allsamples[0]
        sample2 = allsamples[1]
    else:
        if(args.verbose):
            logger.critical("The VCF does not have two sample columns. Please input a proper vcf with Tumor/Normal columns")
        sys.exit(1)
    if(sample1 == args.tsampleName):
        nsampleName = sample2
    else:
        nsampleName = sample1

    # If provided, load hotspots into a dictionary for quick lookup
    hotspot = {}
    if(args.hotspotVcf):
        hvcf_reader = vcf.Reader(open(args.hotspotVcf, 'r'))
        for record in hvcf_reader:
            genomic_locus = str(record.CHROM) + ":" + str(record.POS)
            hotspot[genomic_locus] = True

    for record in vcf_reader:
        tcall = record.genotype(args.tsampleName)
        if(tcall['DP'] is not None):
            tdp = int(tcall['DP'])
        else:
            tdp = 0
        if(tcall['AD'] is not None):
            (trd, tad) = tcall['AD']
        else:
            trd, tad = 0
        if(tdp != 0):
            tvf = int(tad) / float(tdp)
        else:
            tvf = 0
        ncall = record.genotype(nsampleName)
        if(ncall):
            if(ncall['DP'] is not None):
                ndp = int(ncall['DP'])
            else:
                ndp = 0
            if(ncall['AD'] is not None):
                (nrd, nad) = ncall['AD']
            else:
                nad = 0
            if(ndp != 0):
                nvf = nad / ndp
            else:
                nvf = 0
            nvfRF = int(args.tnr) * nvf
        locus = str(record.CHROM) + ":" + str(record.POS)
        if(tvf > nvfRF):
            if((tdp >= int(args.dp)) & (tad >= int(args.ad)) & (tvf >= float(args.vf))):
                vcf_writer.write_record(record)
                txt_fh.write(args.tsampleName + "\t" + record.CHROM + "\t" + str(record.POS) +
                    "\t" + str(record.REF) + "\t" + str(record.ALT[0]) + "\t" + "." + "\n")
        else:
            if(locus in hotspot):
                if((tdp >= int(args.dp)) & (tad >= int(args.ad)) & (tvf >= float(args.vf))):
                    vcf_writer.write_record(record)
                    txt_fh.write(args.tsampleName + "\t" + record.CHROM + "\t" + str(record.POS) +
                        "\t" + str(record.REF) + "\t" + str(record.ALT[0]) + "\t" + "." + "\n")
    vcf_writer.close()
    txt_fh.close()
    return(vcf_out)

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    logging.info("filter_sid: Elapsed time was %g seconds", totaltime)
    sys.exit(0)