#!/usr/bin/env python
'''
@description : This tool helps to filter muTect v1 txt and vcf
@created : 07/17/2016
@author : Ronak H Shah, Cyriac Kandoth

'''

from __future__ import division
import argparse, sys, os, time, logging, cmo, csv

logging.basicConfig(
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logging.DEBUG)
logger = logging.getLogger('filter_mutect')
try:
    import vcf
    from vcf.parser import _Info as VcfInfo, _Format as VcfFormat, _Filter as VcfFilter
except ImportError:
    logger.fatal("filter_mutect: pyvcf is not installed")
    sys.exit(1)

def main():
    parser = argparse.ArgumentParser(prog='filter_mutect.py', description='Filter snps from the output of muTect v1', usage='%(prog)s [options]')
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", help="make lots of noise")
    parser.add_argument("-ivcf", "--inputVcf", action="store", dest="inputVcf", required=True, type=str, metavar='SomeID.vcf', help="Input vcf muTect file which needs to be filtered")
    parser.add_argument("-itxt", "--inputTxt", action="store", dest="inputTxt", required=True, type=str, metavar='SomeID.txt', help="Input txt muTect file which needs to be filtered")
    parser.add_argument("-tsn", "--tsampleName", action="store", dest="tsampleName", required=True, type=str, metavar='SomeName', help="Name of the tumor sample")
    parser.add_argument("-rf", "--refFasta", action="store", dest="refFasta", required=True, type=str, metavar='ref.fa', help="Reference genome in fasta format")
    parser.add_argument("-dp", "--totaldepth", action="store", dest="minDP", required=False, type=int, default=5, metavar='5', help="Minimum Tumor total depth")
    parser.add_argument("-ad", "--alleledepth", action="store", dest="minAD", required=False, type=int, default=3, metavar='3', help="Minimum Tumor allele depth")
    parser.add_argument("-tnr", "--tnRatio", action="store", dest="minTNR", required=False, type=int, default=5, metavar='5', help="Minimum Tumor-Normal variant allele fraction ratio")
    parser.add_argument("-vf", "--variantfraction", action="store", dest="minVAF", required=False, type=float, default=0.01, metavar='0.01', help="Minimum Tumor variant allele fraction")
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
    if(args.outdir):
        vcf_out = os.path.join(args.outdir,vcf_out)
    vcf_out = vcf_out + "_STDfilter.vcf"
    vcf_reader = vcf.Reader(open(args.inputVcf, 'r'))
    vcf_reader.infos['FAILURE_REASON'] = VcfInfo('FAILURE_REASON', '.', 'String', 'Failure Reason from MuTect text File', 'muTect', 'v1')
    vcf_reader.infos['set'] = VcfInfo('set', '.', 'String', 'The variant callers that reported this event', 'mskcc/basicfiltering', 'v0.2.2')
    vcf_reader.formats['DP'] = VcfFormat('DP', '1', 'Integer', 'Total read depth at this site')
    vcf_reader.formats['AD'] = VcfFormat('AD', 'R', 'Integer', 'Allelic depths for the ref and alt alleles in the order listed')
    # Set the soft filter tags we're going to be adding to the VCF
    vaf_tag = 'f' + str(args.minVAF)
    tnr_tag = 'tnr' + str(args.minTNR)
    vcf_reader.filters[vaf_tag] = VcfFilter(vaf_tag, 'Variant Allele Fraction (VAF) <' + str(args.minVAF) + ' in tumor BAM')
    vcf_reader.filters[tnr_tag] = VcfFilter(tnr_tag, 'Non-hotspot with ratio between Tumor-Normal VAFs <' + str(args.minTNR))

    allsamples = list(vcf_reader.samples)
    if len(allsamples) != 2:
        if args.verbose:
            logger.critical("The VCF does not have two genotype columns. Please input a proper vcf with Tumor/Normal columns")
        sys.exit(1)

    # If the caller reported the normal genotype column before the tumor, swap those around
    swap_sample_cols = False
    if allsamples[1] == args.tsampleName:
        swap_sample_cols = True
        vcf_reader.samples[0] = allsamples[1]
        vcf_reader.samples[1] = allsamples[0]
    nsampleName = vcf_reader.samples[1]

    # If provided, load hotspots into a dictionary for quick lookup
    hotspot = {}
    if(args.hotspotVcf):
        hvcf_reader = vcf.Reader(open(args.hotspotVcf, 'r'))
        for record in hvcf_reader:
            genomic_locus = str(record.CHROM) + ":" + str(record.POS)
            hotspot[genomic_locus] = True

    # Parse the MuTect text file to figure out which events to keep
    keepDict = {}
    with open(args.inputTxt, 'rb') as infile:
        reader = csv.DictReader((row for row in infile if not row.startswith('#')), delimiter='\t')
        for row in reader:
            key_for_tracking = str(row['contig']) + ":" + str(row['position']) + ":" + str(row['ref_allele']) + ":" + str(row['alt_allele'])
            rescued_tags = ["alt_allele_in_normal", "nearby_gap_events", "triallelic_site", "possible_contamination", "clustered_read_position"]
            if row['judgement'] == "KEEP" or set(row['failure_reasons'].split(",")).issubset(rescued_tags):
                keepDict[key_for_tracking] = "None" if row['judgement'] == "KEEP" else row['failure_reasons']

    vcf_writer = vcf.Writer(open(vcf_out, 'w'), vcf_reader)
    for record in vcf_reader:
        tcall = record.genotype(args.tsampleName)
        tdp = int(tcall['DP']) if(tcall['DP'] is not None) else 0
        tad = int(tcall['AD'][1]) if(tcall['AD'][1] is not None) else 0
        tvf = int(tad)/float(tdp) if(tdp != 0) else 0
        ncall = record.genotype(nsampleName)
        if ncall:
            ndp = int(ncall['DP']) if(ncall['DP'] is not None) else 0
            nad = int(ncall['AD'][1]) if(ncall['AD'][1] is not None) else 0
            nvf = nad/ndp if(ndp != 0) else 0
            nvfRF = int(args.minTNR) * nvf
        else:
            logger.critical("filter_mutect: There are no genotype values for Normal. We will exit.")
            sys.exit(1)
        locus = str(record.CHROM) + ":" + str(record.POS)
        record.add_info('set', 'MuTect')

        if swap_sample_cols:
            nrm = record.samples[0]
            tum = record.samples[1]
            record.samples[0] = tum
            record.samples[1] = nrm
        key_for_tracking = locus + ":" + str(record.REF) + ":" + str(record.ALT[0])
        if key_for_tracking in keepDict and tdp >= int(args.minDP) and tad >= int(args.minAD):
            # Add some FILTER and INFO tags to the remaining events
            record.FILTER = []
            if tvf < float(args.minVAF):
                record.add_filter(vaf_tag)
            if tvf < nvfRF and locus not in hotspot:
                record.add_filter(tnr_tag)
            record.add_info('FAILURE_REASON', keepDict.get(key_for_tracking))
            vcf_writer.write_record(record)
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
