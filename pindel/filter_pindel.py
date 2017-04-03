#!/usr/bin/python
"""
@Description : This tool helps to filter pindel v0.2.5a7 vcf through command line. 
@Created :  07/17/2014
@Updated: 03/17/2017
@author : Ronak H Shah

"""
from __future__ import division
import argparse
import sys
import time
import os
import logging

logging.basicConfig(
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logging.DEBUG)
logger = logging.getLogger('filter_pindel')
try:
    import coloredlogs
    coloredlogs.install(level='DEBUG')
except ImportError:
    logger.warning("filter_pindel: coloredlogs is not installed, please install it if you wish to see color in logs on standard out.")
    pass
try:
    import vcf
except ImportError:
    logger.fatal("filter_pindel: pyvcf is not installed, please install pyvcf as it is required to run the mapping.")
    sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        prog='filter_pindel.py',
        description='Filter indels from the output of pindel v0.2.5a7',
        usage='%(prog)s [options]')
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        help="make lots of noise")
    parser.add_argument(
        "-ivcf",
        "--inputVcf",
        action="store",
        dest="inputVcf",
        required=True,
        type=str,
        metavar='SomeID.vcf',
        help="Input vcf freebayes file which needs to be filtered")
    parser.add_argument(
        "-tsn",
        "--tsampleName",
        action="store",
        dest="tsampleName",
        required=True,
        type=str,
        metavar='SomeName',
        help="Name of the tumor Sample")
    parser.add_argument(
        "-dp",
        "--totaldepth",
        action="store",
        dest="dp",
        required=False,
        type=int,
        default=0,
        metavar='0',
        help="Tumor total depth threshold")
    parser.add_argument(
        "-ad",
        "--alleledepth",
        action="store",
        dest="ad",
        required=False,
        type=int,
        default=5,
        metavar='5',
        help="Tumor allele depth threshold")
    parser.add_argument(
        "-tnr",
        "--tnRatio",
        action="store",
        dest="tnr",
        required=False,
        type=int,
        default=5,
        metavar='5',
        help="Tumor-Normal variant frequency ratio threshold ")
    parser.add_argument(
        "-vf",
        "--variantfrequency",
        action="store",
        dest="vf",
        required=False,
        type=float,
        default=0.01,
        metavar='0.01',
        help="Tumor variant frequency threshold ")
    parser.add_argument(
        "-o",
        "--outDir",
        action="store",
        dest="outdir",
        required=False,
        type=str,
        metavar='/somepath/output',
        help="Full Path to the output dir.")
    parser.add_argument(
        "-min",
        "--min_var_len",
        action="store",
        dest="min",
        required=False,
        default=25,
        type=int,
        metavar='25',
        help="Minimum length of the indels")
    parser.add_argument(
        "-max",
        "--max_var_len",
        action="store",
        dest="max",
        required=False,
        default=2000,
        type=int,
        metavar='500',
        help="Max length of the indels")
    parser.add_argument(
        "-hvcf",
        "--hotspotVcf",
        action="store",
        dest="hotspotVcf",
        required=False,
        type=str,
        metavar='hostpot.vcf',
        help="Input bgzip / tabix indexed hotspot vcf file to used for filtering")

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
    sample1 = allsamples[0]
    sample2 = allsamples[1]
    if(sample1 == args.tsampleName):
        nsampleName = sample2
    else:
        nsampleName = sample1
    for record in vcf_reader:
        tcall = record.genotype(args.tsampleName)
        recordType = record.INFO['SVTYPE']
        recordLen = abs(int(record.INFO['SVLEN']))
        # print tcall, "recordType: ",recordType, "recordLen: ", recordLen

        if(tcall['AD'] is not None):
            trd, tad = tcall['AD']
        else:
            trd = 0
            tad = 0
        tdp = trd + tad
        if(tdp != 0):
            tvf = int(tad) / float(tdp)
        else:
            tvf = 0
        # print "Tdata: ",trd,tad
        ncall = record.genotype(nsampleName)
        if(ncall):
            if(ncall['AD'] is not None):
                nrd, nad = ncall['AD']
            else:
                nrd = 0
                nad = 0
            ndp = nrd + nad
            if(ndp != 0):
                nvf = nad / ndp
            else:
                nvf = 0
            # print "Ndata: ",trd,tad

            nvfRF = int(args.tnr) * nvf
            # print recordLen, args.min, args.max
        if(args.hotspotVcf):
            hotspotFlag = checkHotspot(args.hotspotVcf, record.CHROM, record.POS)
        else:
            hotspotFlag = False
        if((recordLen >= int(args.min)) and (recordLen <= int(args.max))):
            if(tvf > nvfRF):
                if((tdp >= int(args.dp)) & (tad >= int(args.ad)) & (tvf >= float(args.vf))):
                    if(recordType != "RPL"):
                        # print args.tsampleName , "\t" , record.CHROM , "\t" , str(record.POS) ,
                        # "\t" , str(record.REF) , "\t" + str(record.ALT[0]) , "\t" , "." + "\n"
                        vcf_writer.write_record(record)
                        txt_fh.write(
                            args.tsampleName + "\t" + record.CHROM + "\t" + str(record.POS) + "\t" +
                            str(record.REF) + "\t" + str(record.ALT[0]) + "\t" + "." + "\n")
            else:
                if(hotspotFlag):
                    if((tdp >= int(args.dp)) & (tad >= int(args.ad)) & (tvf >= float(args.vf))):
                        if(recordType != "RPL"):
                            # print args.tsampleName , "\t" , record.CHROM , "\t" , str(record.POS)
                            # , "\t" , str(record.REF) , "\t" + str(record.ALT[0]) , "\t" , "." +
                            # "\n"
                            vcf_writer.write_record(record)
                            txt_fh.write(args.tsampleName + "\t" + record.CHROM + "\t" +
                                         str(record.POS) + "\t" + str(record.REF) + "\t" +
                                         str(record.ALT[0]) + "\t" + "." + "\n")

    txt_fh.close()
    return(vcf_out)

def checkHotspot(hotspotVcf, chromosome, start):
    hotspotFlag = False
    hotspot_vcf_reader = vcf.Reader(open(hotspotVcf, 'r'))
    try:
        record = hotspot_vcf_reader.fetch(str(chromosome), start)
    except ValueError:
        logger.info("filter_pindel: Region not present in vcf, %s:%s", str(chromosome), start)
        record = None

    if(record is None):
        hotspotFlag = False
    else:
        hotspotFlag = True
    return(hotspotFlag)

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    logging.info("filter_pindel: Elapsed time was %g seconds", totaltime)
    sys.exit(0)
