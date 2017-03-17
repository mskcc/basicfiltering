'''
@Description : This tool helps to filter SomaticIndelDetector vcf and txt through command line.
@Created :  07/25/2016
@author : Ronak H Shah

'''
from __future__ import division
import argparse
import sys
import os
import time
import os.path
import vcf


def main():
    parser = argparse.ArgumentParser(
        prog='filter_sid.py',
        description='Filter Indels from the output of SomaticIndelDetector',
        usage='%(prog)s [options]')
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        default=True,
        help="make lots of noise [default]")
    parser.add_argument(
        "-ivcf",
        "-inputVcf",
        action="store",
        dest="inputVcf",
        required=True,
        metavar='SomeID.vcf',
        help="Input SomaticIndelDetector vcf file which needs to be filtered")
    parser.add_argument(
        "-itxt",
        "-inputTxt",
        action="store",
        dest="inputTxt",
        required=True,
        metavar='SomeID.txt',
        help="Input SomaticIndelDetector txt file which needs to be filtered")
    parser.add_argument(
        "-tsn",
        "--tsampleName",
        action="store",
        dest="tsampleName",
        required=True,
        metavar='SomeName',
        help="Name of the tumor Sample")
    parser.add_argument(
        "-dp",
        "--totaldepth",
        action="store",
        dest="dp",
        required=False,
        default=0,
        metavar='0',
        help="Tumor total depth threshold")
    parser.add_argument(
        "-ad",
        "--alleledepth",
        action="store",
        dest="ad",
        required=False,
        default=5,
        metavar='5',
        help="Tumor allele depth threshold")
    parser.add_argument(
        "-tnr",
        "--tnRatio",
        action="store",
        dest="tnr",
        required=False,
        default=5,
        metavar='5',
        help="Tumor-Normal variant frequency ratio threshold ")
    parser.add_argument(
        "-vf",
        "--variantfrequency",
        action="store",
        dest="vf",
        required=False,
        default=0.01,
        metavar='0.01',
        help="Tumor variant frequency threshold ")
    parser.add_argument(
        "-hvcf",
        "--hotspotVcf",
        action="store",
        dest="hotspotVcf",
        required=True,
        metavar='hostpot.vcf',
        help="Input bgzip / tabix indexed hotspot vcf file to used for filtering")
    parser.add_argument(
        "-o",
        "--outDir",
        action="store",
        dest="outdir",
        required=True,
        metavar='/somepath/output',
        help="Full Path to the output dir.")

    args = parser.parse_args()
    if(args.verbose):
        print "I have Started the run for doing standard filter."
    (stdfilterVCF) = RunStdFilter(args)
    if(args.verbose):
        print "I have finished the run for doing standard filter."


# Code that does Standard Filter
def RunStdFilter(args):
    vcf_out = os.path.basename(args.inputVcf)
    vcf_out = os.path.splitext(vcf_out)[0]
    txt_out = vcf_out
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
        # print "Tdata: ",trd,tad
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
            # print "Ndata: ",trd,tad

            nvfRF = int(args.tnr) * nvf
        hotspotFlag = checkHotspot(args.hotspotVcf, record.CHROM, record.POS)
        if(tvf > nvfRF):
            if((tdp >= int(args.dp)) & (tad >= int(args.ad)) & (tvf >= float(args.vf))):
                vcf_writer.write_record(record)
                txt_fh.write(args.tsampleName +
                             "\t" +
                             record.CHROM +
                             "\t" +
                             str(record.POS) +
                             "\t" +
                             str(record.REF) +
                             "\t" +
                             str(record.ALT[0]) +
                             "\t" +
                             "." +
                             "\n")

        else:
            if(hotspotFlag):
                if((tdp >= int(args.dp)) & (tad >= int(args.ad)) & (tvf >= float(args.vf))):
                    vcf_writer.write_record(record)
                    txt_fh.write(
                        args.tsampleName + "\t" + record.CHROM + "\t" + str(record.POS) + "\t" +
                        str(record.REF) + "\t" + str(record.ALT[0]) + "\t" + "." + "\n")

    txt_fh.close()
    return(vcf_out)


def checkHotspot(hotspotVcf, chromosome, start):
    hotspotFlag = False
    hotspot_vcf_reader = vcf.Reader(open(hotspotVcf, 'r'))
    try:
        record = hotspot_vcf_reader.fetch(str(chromosome), start)
    except ValueError:
        print ("Region not present in vcf, ", str(chromosome), ":", start)
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
    print("Elapsed time was %g seconds" % (end_time - start_time))
