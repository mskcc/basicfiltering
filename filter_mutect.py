#!/usr/bin/env python
'''
@description : This tool helps to filter muTect v1 txt and vcf
@created : 07/17/2016
@author : Ronak H Shah, Cyriac Kandoth

'''

from __future__ import division
import argparse, sys, os, time, logging, cmo

logging.basicConfig(
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logging.DEBUG)
logger = logging.getLogger('filter_mutect')
try:
    import vcf
    from vcf.parser import _Info as VcfInfo, _Format as VcfFormat
except ImportError:
    logger.fatal("filter_mutect: pyvcf is not installed")
    sys.exit(1)
try:
    import pandas as pd
except ImportError:
    logger.fatal("filter_mutect: pandas is not installed")
    sys.exit(1)

def main():
    parser = argparse.ArgumentParser(prog='filter_mutect.py', description='Filter snps from the output of muTect v1', usage='%(prog)s [options]')
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", help="make lots of noise")
    parser.add_argument("-ivcf", "--inputVcf", action="store", dest="inputVcf", required=True, type=str, metavar='SomeID.vcf', help="Input vcf muTect file which needs to be filtered")
    parser.add_argument("-itxt", "--inputTxt", action="store", dest="inputTxt", required=True, type=str, metavar='SomeID.txt', help="Input txt muTect file which needs to be filtered")
    parser.add_argument("-tsn", "--tsampleName", action="store", dest="tsampleName", required=True, type=str, metavar='SomeName', help="Name of the tumor Sample")
    parser.add_argument("-rf", "--refFasta", action="store", dest="refFasta", required=True, type=str, metavar='ref.fa', help="Reference genome in fasta format")
    parser.add_argument("-dp", "--totaldepth", action="store", dest="dp", required=False, type=int, default=5, metavar='5', help="Tumor total depth threshold")
    parser.add_argument("-ad", "--alleledepth", action="store", dest="ad", required=False, type=int, default=3, metavar='3', help="Tumor allele depth threshold")
    parser.add_argument("-tnr", "--tnRatio", action="store", dest="tnr", required=False, type=int, default=5, metavar='5', help="Tumor-Normal variant fraction ratio threshold ")
    parser.add_argument("-vf", "--variantfraction", action="store", dest="vf", required=False, type=float, default=0.01, metavar='0.01', help="Tumor variant fraction threshold ")
    parser.add_argument("-hvcf", "--hotspotVcf", action="store", dest="hotspotVcf", required=False, type=str, metavar='hotspot.vcf', help="Input vcf file with hotspots that skip VAF ratio filter")
    parser.add_argument("-o", "--outDir", action="store", dest="outdir", required=False, type=str, metavar='/somepath/output', help="Full Path to the output dir.")

    args = parser.parse_args()
    if(args.verbose):
        logger.info("Started the run for doing standard filter.")
    RunStdFilter(args)
    if(args.verbose):
        logger.info("Finished the run for doing standard filter.")

def RunStdFilter(args):
    vcf_out = os.path.basename(args.inputVcf)
    vcf_out = os.path.splitext(vcf_out)[0]
    txt_out = os.path.basename(args.inputTxt)
    txt_out = os.path.splitext(txt_out)[0]
    if(args.outdir):
        vcf_out = os.path.join(args.outdir,vcf_out)
        txt_out = os.path.join(args.outdir,txt_out)
    vcf_out = vcf_out + "_STDfilter.vcf"
    txt_out = txt_out + "_STDfilter.txt"
    vcf_reader = vcf.Reader(open(args.inputVcf, 'r'))
    vcf_reader.infos['FAILURE_REASON'] = VcfInfo('FAILURE_REASON', '.', 'String', 'Failure Reason from MuTect text File', 'muTect', 'v1')
    vcf_reader.infos['set'] = VcfInfo('set', '.', 'String', 'The variant callers that reported this event', 'mskcc/basicfiltering', 'v0.2.2')
    vcf_reader.formats['DP'] = VcfFormat('DP', '1', 'Integer', 'Total read depth at this site')
    vcf_reader.formats['AD'] = VcfFormat('AD', 'R', 'Integer', 'Allelic depths for the ref and alt alleles in the order listed')

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

    # Dictionary to store records to keep
    keepDict = {}

    # If provided, load hotspots into a dictionary for quick lookup
    hotspot = {}
    if(args.hotspotVcf):
        hvcf_reader = vcf.Reader(open(args.hotspotVcf, 'r'))
        for record in hvcf_reader:
            genomic_locus = str(record.CHROM) + ":" + str(record.POS)
            hotspot[genomic_locus] = True

    txtDF = pd.read_table(args.inputTxt, skiprows=1, dtype=str)
    txt_fh = open(txt_out, "wb")
    for index, row in txtDF.iterrows():
        chr = row.loc['contig']  # Get Chromosome
        pos = row.loc['position']  # Get Position
        ref_allele = row.loc['ref_allele']
        alt_allele = row.loc['alt_allele']
        trd = int(row.loc['t_ref_count'])
        tad = int(row.loc['t_alt_count'])
        tdp = trd + tad
        tvf = int(tad)/float(tdp) if(tdp != 0) else 0
        nrd = int(row.loc['n_ref_count'])
        nad = int(row.loc['n_alt_count'])
        ndp = nrd + nad
        nvf = int(nad) / float(ndp) if(ndp != 0) else 0
        judgement = row.loc['judgement']  # Get REJECT or PASS
        failure_reason = row.loc['failure_reasons']  # Get Reject Reason
        nvfRF = int(args.tnr) * nvf

        # This will help in filtering VCF
        key_for_tracking = str(chr) + ":" + str(pos) + ":" + str(ref_allele) + ":" + str(alt_allele)
        locus = str(chr) + ":" + str(pos)
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
            if tvf > nvfRF or locus in hotspot:
                if((tdp >= int(args.dp)) & (tad >= int(args.ad)) & (tvf >= float(args.vf))):
                    if(key_for_tracking in keepDict):
                        print("MutectStdFilter:There is a repeat ", key_for_tracking)
                    else:
                        keepDict[key_for_tracking] = failure_reason
                    txt_fh.write(args.tsampleName + "\t" + str(chr) + "\t" + str(pos) + "\t" + str(ref_allele) + "\t" + str(alt_allele) + "\t" + str(failure_reason) + "\n")
    txt_fh.close()

    vcf_writer = vcf.Writer(open(vcf_out, 'w'), vcf_reader)
    for record in vcf_reader:
        key_for_tracking = str(record.CHROM) + ":" + str(record.POS) + ":" + str(record.REF) + ":" + str(record.ALT[0])
        if(key_for_tracking in keepDict):
            failure_reason = keepDict.get(key_for_tracking)
            if(failure_reason == "KEEP"):
                failure_reason = "None"
            record.add_info('FAILURE_REASON', failure_reason)
            record.add_info('set', 'MuTect')
            if (if_swap_sample):
                nrm = record.samples[0]
                tum = record.samples[1]
                record.samples[0] = tum
                record.samples[1] = nrm
            if(record.FILTER == "PASS"):
                vcf_writer.write_record(record)
            else:
                record.FILTER = "PASS"
                vcf_writer.write_record(record)

        else:
            continue
    vcf_writer.close()
    # Normalize the events in the VCF, produce a bgzipped VCF, then tabix index it
    norm_gz_vcf = cmo.util.normalize_vcf(vcf_out, args.refFasta)
    cmo.util.tabix_file(norm_gz_vcf)
    return(norm_gz_vcf)

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    logging.info("filter_mutect: Elapsed time was %g seconds", totaltime)
    sys.exit(0)
