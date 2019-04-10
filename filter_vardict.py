#!/usr/bin/env python
'''
@description : Filter a VCF listing somatic events called by VarDict
@created : 04/22/2016
@author : Ronak H Shah, Cyriac Kandoth, Zuojian Tang

'''

from __future__ import division
import argparse, sys, os, time, logging, cmo

logging.basicConfig(
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logging.DEBUG)
logger = logging.getLogger('filter_vardict')
try:
    import vcf
    from vcf.parser import _Info as VcfInfo, _Format as VcfFormat, _Filter as VcfFilter
except ImportError:
    logger.fatal("filter_mutect: pyvcf is not installed")
    sys.exit(1)

def main():
    parser = argparse.ArgumentParser(prog='filter_vardict.py',description='Filter a VCF listing somatic events called by VarDict',usage='%(prog)s [options]')
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",help="make lots of noise")
    parser.add_argument("-ivcf", "--inputVcf", action="store", dest="inputVcf", required=True, type=str, metavar='SomeID.vcf',help="Input vcf vardict file which needs to be filtered")
    parser.add_argument("-tsn", "--tsampleName", action="store", dest="tsampleName", required=True, type=str, metavar='SomeName',help="Name of the tumor Sample")
    parser.add_argument("-rf", "--refFasta", action="store", dest="refFasta", required=True, type=str, metavar='ref.fa', help="Reference genome in fasta format")
    parser.add_argument("-dp", "--totaldepth", action="store", dest="minDP", required=False, type=int, default=5, metavar='5',help="Tumor total depth threshold")
    parser.add_argument("-ad", "--alleledepth", action="store", dest="minAD", required=False, type=int, default=3, metavar='3',help="Tumor allele depth threshold")
    parser.add_argument("-tnr", "--tnRatio", action="store", dest="minTNR", required=False, type=int, default=5, metavar='5',help="Tumor-Normal variant fraction ratio threshold ")
    parser.add_argument("-vf", "--variantfraction", action="store", dest="minVAF", required=False, type=float, default=0.01, metavar='0.01',help="Tumor variant fraction threshold ")
    parser.add_argument("-mq", "--minqual", action="store", dest="minQual", required=False, type=int, default=20, metavar='20',help="Minimum variant call quality")
    parser.add_argument("-cm", "--cpxMinMQ", action="store", dest="cpxMinMQ", required=False, type=float, default=55, metavar='55',help="Minimum mean mapping qual for complex events")
    parser.add_argument("-cn", "--cpxMaxNM", action="store", dest="cpxMaxNM", required=False, type=float, default=2, metavar='2',help="Maximum mean mismatches for complex events")
    parser.add_argument("-hvcf", "--hotspotVcf", action="store", dest="hotspotVcf", required=False, type=str, metavar='hotspot.vcf',help="Input vcf file with hotspots that skip VAF ratio filter")
    parser.add_argument("-o", "--outDir", action="store", dest="outdir", required=False, type=str, metavar='/somepath/output',help="Full Path to the output dir.")

    args = parser.parse_args()
    if args.verbose:
        logger.info("Started the run for doing standard filter.")
    RunStdFilter(args)
    if args.verbose:
        logger.info("Finished the run for doing standard filter.")

def RunStdFilter(args):
    vcf_out = os.path.basename(args.inputVcf)
    vcf_out = os.path.splitext(vcf_out)[0]
    if args.outdir:
        vcf_out = os.path.join(args.outdir,vcf_out)
    vcf_out = vcf_out + "_STDfilter.vcf"
    vcf_reader = vcf.Reader(open(args.inputVcf, 'r'))
    vcf_reader.infos['set'] = VcfInfo('set', '.', 'String', 'The variant callers that reported this event', 'mskcc/basicfiltering', 'v0.2.2')
    vcf_reader.infos['VSB'] = VcfInfo('VSB', '0', 'Flag', 'Non-hotspot with strand-bias (VD>10 and ALD lists 0) in tumor data, unless all REF/ALT reads have strand-bias at MQ>40', 'mskcc/basicfiltering', 'v0.2.2')
    vcf_reader.formats['DP'] = VcfFormat('DP', '1', 'Integer', 'Total read depth at this site')
    vcf_reader.formats['AD'] = VcfFormat('AD', 'R', 'Integer', 'Allelic depths for the ref and alt alleles in the order listed')
    # Set the soft filter tags we're going to be adding to the VCF
    vaf_tag = 'f' + str(args.minVAF)
    tnr_tag = 'tnr' + str(args.minTNR)
    nm_tag = 'nm' + str(args.cpxMaxNM)
    mq_tag = 'mq' + str(args.cpxMinMQ)
    vcf_reader.filters[vaf_tag] = VcfFilter(vaf_tag, 'Variant Allele Fraction (VAF) <' + str(args.minVAF) + ' in tumor BAM')
    vcf_reader.filters[tnr_tag] = VcfFilter(tnr_tag, 'Non-hotspot with ratio between Tumor-Normal VAFs <' + str(args.minTNR))
    vcf_reader.filters[nm_tag] = VcfFilter(nm_tag, 'Non-hotspot non-SNV with a mean number of mismatches >' + str(args.cpxMaxNM) + ' in tumor BAM')
    vcf_reader.filters[mq_tag] = VcfFilter(mq_tag, 'Non-hotspot non-SNV with a mean mapping quality <' + str(args.cpxMinMQ) + ' in tumor BAM')
    # Set hstdp, hsndp, hstad, hsvaf tags
    hstdp_tag = 'hstdp'
    hsndp_tag = 'hsndp'
    hstad_tag = 'hstad'
    hsvaf_tag = 'hsvaf'
    vcf_reader.filters[hstdp_tag] = VcfFilter(hstdp_tag, 'Tumor depth <12 for hotspots, or <20 for non-hotspots')
    vcf_reader.filters[hsndp_tag] = VcfFilter(hsndp_tag, 'Normal depth <6 for hotspots, or <10 for non-hotspots')
    vcf_reader.filters[hstad_tag] = VcfFilter(hstad_tag, 'Tumor allele depth <3 for hotspots, or <5 for non-hotspots')
    vcf_reader.filters[hsvaf_tag] = VcfFilter(hsvaf_tag, 'Variant allele fraction <0.02 for hotspots, or <0.05 for non-hotspots')

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
    if args.hotspotVcf:
        hvcf_reader = vcf.Reader(open(args.hotspotVcf, 'r'))
        for record in hvcf_reader:
            genomic_locus = str(record.CHROM) + ":" + str(record.POS)
            hotspot[genomic_locus] = True

    vcf_writer = vcf.Writer(open(vcf_out, 'w'), vcf_reader)
    for record in vcf_reader:
        tcall = record.genotype(args.tsampleName)
        somaticStatus = True if "Somatic" in record.INFO['STATUS'] else False
        tql = int(tcall['QUAL']) if(tcall['QUAL'] is not None) else 0
        tdp = int(tcall['DP']) if(tcall['DP'] is not None) else 0
        tad = int(tcall['VD']) if(tcall['VD'] is not None) else 0
        tnm = float(tcall['NM']) if(tcall['NM'] is not None) else 0
        tmq = float(tcall['MQ']) if(tcall['MQ'] is not None) else 0
        tvf = int(tad)/float(tdp) if(tdp != 0) else 0
        if tcall['ALD'] and tcall['RD']:
            tdp_fwdrev = [tcall['ALD'][0] + tcall['RD'][0], tcall['ALD'][1] + tcall['RD'][1]]

        ncall = record.genotype(nsampleName)
        if ncall:
            nql = int(ncall['QUAL']) if(ncall['QUAL'] is not None) else 0
            ndp = int(ncall['DP']) if(ncall['DP'] is not None) else 0
            nad = int(ncall['VD']) if(ncall['VD'] is not None) else 0
            nvf = nad/ndp if(ndp != 0) else 0
            nvfRF = int(args.minTNR) * nvf
        else:
            logger.critical("filter_vardict: There are no genotype values for Normal. We will exit.")
            sys.exit(1)
        locus = str(record.CHROM) + ":" + str(record.POS)
        record.add_info('set', 'VarDict')

        if swap_sample_cols:
            nrm = record.samples[0]
            tum = record.samples[1]
            record.samples[0] = tum
            record.samples[1] = nrm
        if somaticStatus and tql >= int(args.minQual) and nql >= int(args.minQual) and tdp >= int(args.minDP) and tad >= int(args.minAD):
            # Add some FILTER and INFO tags to the remaining events
            if tvf < float(args.minVAF):
                record.add_filter(vaf_tag)
            if tvf < nvfRF and locus not in hotspot:
                record.add_filter(tnr_tag)
            if locus not in hotspot and record.INFO['TYPE'] != "SNV":
                if tnm > args.cpxMaxNM:
                    record.add_filter(nm_tag)
                if tmq < args.cpxMinMQ:
                    record.add_filter(mq_tag)
            # Non-hotspot with strand-bias (VD>10 and ALD lists 0) in tumor data, unless all REF/ALT reads have strand-bias at MQ>40
            if locus not in hotspot and tad >= 10 and 0 in tcall['ALD'] and not (0 in tdp_fwdrev and tmq > 40):
                record.add_info('VSB')
            if (tdp < 12 and locus in hotspot) or (tdp < 20 and locus not in hotspot):
                record.add_filter(hstdp_tag)
            if (ndp < 6 and locus in hotspot) or (ndp < 10 and locus not in hotspot):
                record.add_filter(hsndp_tag)
            if (tad < 3 and locus in hotspot) or (tad < 5 and locus not in hotspot):
                record.add_filter(hstad_tag)
            if (tvf < 0.02 and locus in hotspot) or (tvf < 0.05 and locus not in hotspot):
                record.add_filter(hsvaf_tag)
            vcf_writer.write_record(record)
    vcf_writer.close()
    # Normalize the events in the VCF, produce a bgzipped VCF, then tabix index it
    norm_gz_vcf = cmo.util.normalize_vcf(vcf_out, args.refFasta)
    cmo.util.tabix_file(norm_gz_vcf)
    return norm_gz_vcf

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    logging.info("filter_vardict: Elapsed time was %g seconds", totaltime)
    sys.exit(0)
