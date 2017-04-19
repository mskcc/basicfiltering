#!/usr/bin/python
'''
@Description : This tool helps to filter muTect v1.14 txt and vcf through command line. 
@Created :  07/17/2016
@Updated : 03/17/2017
@author : Ronak H Shah

'''
from __future__ import division
import argparse
import sys
import os
import time
import logging

logging.basicConfig(
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logging.DEBUG)
logger = logging.getLogger('filter_mutect')
try:
    import coloredlogs
    coloredlogs.install(level='DEBUG')
except ImportError:
    logger.warning("filter_mutect: coloredlogs is not installed, please install it if you wish to see color in logs on standard out.")
    pass
try:
    import vcf
    from vcf.parser import _Info as VcfInfo, field_counts as vcf_field_counts
except ImportError:
    logger.fatal("filter_mutect: pyvcf is not installed, please install pyvcf as it is required to run the mapping.")
    sys.exit(1)
try:
    import pandas as pd
except ImportError:
    logger.fatal("filter_mutect: pandas is not installed, please install pandas as it is required to run the mapping.")
    sys.exit(1)

def main():
   parser = argparse.ArgumentParser(prog='filter_mutect.py', description='Filter snps from the output of muTect v1.14', usage='%(prog)s [options]')
   parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", help="make lots of noise")
   parser.add_argument("-ivcf", "--inputVcf", action="store", dest="inputVcf", required=True, type=str, metavar='SomeID.vcf', help="Input vcf muTect file which needs to be filtered")
   parser.add_argument("-itxt", "--inputTxt", action="store", dest="inputTxt", required=True, type=str, metavar='SomeID.txt', help="Input txt muTect file which needs to be filtered")
   parser.add_argument("-tsn", "--tsampleName", action="store", dest="tsampleName", required=True, type=str,metavar='SomeName', help="Name of the tumor Sample")
   parser.add_argument("-dp", "--totaldepth", action="store", dest="dp", required=False, type=int, default=0, metavar='0', help="Tumor total depth threshold")
   parser.add_argument("-ad", "--alleledepth", action="store", dest="ad", required=False, type=int, default=5, metavar='5', help="Tumor allele depth threshold")
   parser.add_argument("-tnr", "--tnRatio", action="store", dest="tnr", required=False, type=int, default=5, metavar='5', help="Tumor-Normal variant frequency ratio threshold ")
   parser.add_argument("-vf", "--variantfrequency", action="store", dest="vf", required=False, type=float, default=0.01, metavar='0.01', help="Tumor variant frequency threshold ")
   parser.add_argument("-hvcf", "--hotspotVcf", action="store", dest="hotspotVcf", required=False, type=str, metavar='hostpot.vcf', help="Input bgzip / tabix indexed hotspot vcf file to used for filtering")
   parser.add_argument("-o", "--outDir", action="store", dest="outdir", required=False, type=str, metavar='/somepath/output', help="Full Path to the output dir.")
   
   args = parser.parse_args()
   if(args.verbose):
       logger.info("Started the run for doing standard filter.")
   (stdfilterVCF) = RunStdFilter(args)
   if(args.verbose):
       logger.info("Finished the run for doing standard filter.")
       
# Code that does Standard Filter  
def RunStdFilter(args):
    vcf_out = os.path.basename(args.inputVcf)
    vcf_out = os.path.splitext(vcf_out)[0]
    txt_out = os.path.basename(args.inputTxt)
    txt_out = os.path.splitext(txt_out)[0]
    if(args.outdir):
        vcf_out = os.path.join(args.outdir,vcf_out + "_STDfilter.vcf")
        txt_out = os.path.join(args.outdir,txt_out + "_STDfilter.txt")
    else:
        vcf_out = vcf_out + "_STDfilter.vcf"
        txt_out = txt_out + "_STDfilter.txt"
    vcf_reader = vcf.Reader(open(args.inputVcf, 'r'))   
    vcf_reader.infos['FAILURE_REASON'] = VcfInfo('FAILURE_REASON', 'String', "1", 'Failure Reason from MuTect text File')
    vcf_writer = vcf.Writer(open(vcf_out, 'w'), vcf_reader)
    txtDF = pd.read_table(args.inputTxt, skiprows=1, low_memory=False)
    txt_fh = open(txt_out, "wb") 
    allsamples = vcf_reader.samples
    sample1 = allsamples[0]
    sample2 = allsamples[1]
    if(sample1 == args.tsampleName):
        nsampleName = sample2
    else:
        nsampleName = sample1 
    
    # Dictionalry to store records to keep
    keepDict = {}

    for index, row in txtDF.iterrows():
        chr = row.loc['contig']  # Get Chromosome
        pos = row.loc['position']  # Get Position
        ref_allele = row.loc['ref_allele'] 
        alt_allele = row.loc['alt_allele']
        trd = int(row.loc['t_ref_count'])
        tad = int(row.loc['t_alt_count'])
        tdp = trd + tad 
        if(tdp != 0):
            tvf = int(tad) / float(tdp)
        else:
            tvf = 0
        nrd = int(row.loc['n_ref_count'])
        nad = int(row.loc['n_alt_count'])
        ndp = nrd + nad 
        if(ndp != 0):
            nvf = int(nad) / float(ndp)
        else:
            nvf = 0
        judgement = row.loc['judgement']  # Get REJECT or PASS
        failure_reason = row.loc['failure_reasons']  # Get Reject Reason
        nvfRF = int(args.tnr) * nvf 
        if(args.hotspotVcf):
            hotspotFlag = checkHotspot(args.hotspotVcf, chr, pos) 
        else:
            hotspotFlag = False
        
        # This will help in filtering VCF
        key_for_tracking = str(chr) + ":" + str(pos) + ":" + str(ref_allele) + ":" + str(alt_allele)
        if(judgement == "KEEP"):
            
            if(key_for_tracking in keepDict):
                print("MutectStdFilter:There is a repeat ", key_for_tracking)
            else:
                keepDict[key_for_tracking] = judgement
                txt_fh.write(args.tsampleName + "\t" + str(chr) + "\t" + str(pos) + "\t" + str(ref_allele) + "\t" + str(alt_allele) + "\t" + str(judgement) + "\n")
        else:
            accepted_tags = ["alt_allele_in_normal", "nearby_gap_events", "triallelic_site", "possible_contamination", "clustered_read_position"]
            failure_tags = failure_reason.split(",")
            tag_count = 0
            for tag in failure_tags:
                if tag in accepted_tags:
                    tag_count = tag_count + 1
                else:
                    continue
            if(tag_count != len(failure_tags)):
                continue
            if(tvf > nvfRF):
                if((tdp >= int(args.dp)) & (tad >= int(args.ad)) & (tvf >= float(args.vf))):
                    if(key_for_tracking in keepDict):
                        print("MutectStdFilter:There is a repeat ", key_for_tracking)
                    else:
                        keepDict[key_for_tracking] = failure_reason
                    txt_fh.write(args.tsampleName + "\t" + str(chr) + "\t" + str(pos) + "\t" + str(ref_allele) + "\t" + str(alt_allele) + "\t" + str(failure_reason) + "\n")
            else:
                if(hotspotFlag):
                    if((tdp >= int(args.dp)) & (tad >= int(args.ad)) & (tvf >= float(args.vf))):
                        if(key_for_tracking in keepDict):
                            print("MutectStdFilter:There is a repeat ", key_for_tracking)
                        else:
                            keepDict[key_for_tracking] = failure_reason
                        txt_fh.write(args.tsampleName + "\t" + str(chr) + "\t" + str(pos) + "\t" + str(ref_allele) + "\t" + str(alt_allele) + "\t" + str(failure_reason) + "\n")
    
    txt_fh.close()
    for record in vcf_reader:
        key_for_tracking = str(record.CHROM) + ":" + str(record.POS) + ":" + str(record.REF) + ":" + str(record.ALT[0]) 
        if(key_for_tracking in keepDict):
            failure_reason = keepDict.get(key_for_tracking)
            record.add_info('FAILURE_REASON', failure_reason)
            if(record.FILTER == "PASS"):
                 
                vcf_writer.write_record(record)
            else:
                record.FILTER = "PASS"
                vcf_writer.write_record(record)
               
        else:
            continue
    vcf_writer.close()
    return(vcf_out)

def checkHotspot(hotspotVcf, chromosome, start):
    hotspotFlag = False
    hotspot_vcf_reader = vcf.Reader(open(hotspotVcf, 'r'))
    try:
        record = hotspot_vcf_reader.fetch(str(chromosome), start)
    except ValueError:
        logger.info("filter_mutect: Region not present in vcf, %s:%s", str(chromosome), start)
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
    logging.info("filter_mutect: Elapsed time was %g seconds", totaltime)
    sys.exit(0)
