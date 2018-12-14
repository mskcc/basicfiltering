#!/usr/bin/python

'''
@Description : This tool helps to post filtering Vardict's complex/Insertion/Deletion events by retrieving indels/soft-clipping information from T/N BAM files
@Created :  11/12/2018
'''

from __future__ import division
from __future__ import print_function

import argparse, time, os, sys, logging, re, csv, glob, subprocess
import pysam
from pysam import VariantFile
import inspect

def main():

    # some parameters specific to VarDict output
    str_cplx = "Complex"
    str_dele = "Deletion"
    str_isrt = "Insertion"

    parser = argparse.ArgumentParser(prog='basic-filtering.vardict-complex-events.py',
                                     description='This tool helps to post filtering Vardict\'s complex events by retrieving indels/soft-clipping information from T/N BAM files.',
                                     usage='%(prog)s [options]')
    parser.add_argument("-VCF", "--vcffile", action="store", dest="vcffile", required=True, type=str,
                        help="VarDict VCF file")
    parser.add_argument("-T", "--tumorbam", action="store", dest="tumorbam", required=True, type=str,
                        help="tumor bam file")
    parser.add_argument("-N", "--normalbam", action="store", dest="normalbam", required=True, type=str,
                        help="normal bam file")
    parser.add_argument('-TN', '--tumorname', action="store", dest="tumorname", required=True, type=str,
                        help="tumor sample ID")
    parser.add_argument("-o", "--output", action="store", dest="output", required=True, type=str,
                        help="Output file name including path")
    parser.add_argument("-LEN", "--flanklen", action="store", dest="flanklen", required=False,
                        type=int, default=50, help="flanking length around position [default: 50]")
    parser.add_argument("-MQ", "--mappingquality", action="store", dest="mappingquality", required=False,
                        type=int, default=20, help="mapping quality [default: 20]")
    parser.add_argument("-TNOISE", "--tnoise", action="store", dest="tnoise", required=False,
                        type=float, default=0.2, help="tumor noise (number of reads with indels or soft-clipping/total depth reported by VarDict [default: 0.2]")
    parser.add_argument("-NNoise", "--nnoise", action="store", dest="nnoise", required=False,
                        type=float, default=0.1, help="normal noise (number of reads with indels or soft-clipping/total depth reported by VarDict [default: 0.1]")

    args = parser.parse_args()
    vcf_in = args.vcffile
    t = args.tumorbam
    n = args.normalbam
    tn = args.tumorname
    flen = args.flanklen
    mqual = args.mappingquality
    tnoise = args.tnoise
    nnoise = args.nnoise

    # figure output vcf file
    vcf_out = args.output
    vcf_out_dir = os.path.abspath(vcf_out)
    vcf_out_forR = vcf_out_dir + ".score.forR"
    outr = vcf_out_forR + ".R"
    out_forR = []

    # read bam file
    tbam = pysam.AlignmentFile(t, "rb")
    nbam = pysam.AlignmentFile(n, "rb")

    # read vcf file
    vcf_in_fr = VariantFile(vcf_in)
    vcf_in_fr.header.formats.add("IS", "2", "Integer",
                                 "(Number of reads with indels, number of reads with soft-clipping) [within the flanking region of event]")
    vcf_in_fr.header.filters.add("CPX", None, None, "VarDict Complex and Indels events.")
    vcf_out_fw = VariantFile(vcf_out, 'w', header=vcf_in_fr.header)

    allsamples = vcf_in_fr.header.samples
    if len(allsamples) == 2:
        tid = allsamples[0]
        nid = allsamples[1]
    else:
        print("The VCF does not have two sample columns. Please input a proper vcf with Tumor/Normal columns")
        sys.exit(1)
    if_swap_sample = False
    if nid == tn:
        if_swap_sample = True

    for vcf_in_row in vcf_in_fr.fetch():
        chr = vcf_in_row.chrom
        pos = vcf_in_row.pos
        lpos = pos - flen
        if lpos < 0:
            lpos = 0
        rpos = pos + flen
        ref = vcf_in_row.ref
        alts = vcf_in_row.alts
        tid_tmp = vcf_in_row.samples.keys()[0]
        nid_tmp = vcf_in_row.samples.keys()[1]
        tcall = vcf_in_row.samples[0]
        ncall = vcf_in_row.samples[1]
        if if_swap_sample:
            tid_tmp = vcf_in_row.samples.keys()[1]
            nid_tmp = vcf_in_row.samples.keys()[0]
            tcall = vcf_in_row.samples[1]
            ncall = vcf_in_row.samples[0]
        if tcall['DP'] is not None:
            tdp = tcall['DP']
        if tcall['VD'] is not None:
            tvd = tcall['VD']
        if ncall['DP'] is not None:
            ndp = ncall['DP']
        if ncall['VD'] is not None:
            nvd = ncall['VD']
        type = vcf_in_row.info["TYPE"]
        tcounter_reads_sf = 0
        ncounter_reads_sf = 0
        tcounter_reads_indels = 0
        ncounter_reads_indels = 0
        for tread in tbam.fetch(chr, lpos, rpos):
            tifsrd = tread.is_secondary
            if tifsrd:
                continue
            tifdup = tread.is_duplicate
            if tifdup:
                continue
            tmqual = tread.mapping_quality
            if tmqual < mqual:
                continue
            tcigstr = tread.cigarstring
            if tcigstr is None:
                continue
            if 'I' in tcigstr or 'D' in tcigstr:
                tcounter_reads_indels += 1
            elif 'S' in tcigstr:
                tcounter_reads_sf += 1
        for nread in nbam.fetch(chr, lpos, rpos):
            nifsrd = nread.is_secondary
            if nifsrd:
                continue
            nifdup = nread.is_duplicate
            if nifdup:
                continue
            nmqual = nread.mapping_quality
            if nmqual < mqual:
                continue
            ncigstr = nread.cigarstring
            if ncigstr is None:
                continue
            if 'I' in ncigstr or 'D' in ncigstr:
                ncounter_reads_indels += 1
            elif 'S' in ncigstr:
                ncounter_reads_sf += 1
        vcf_in_row.samples[tid_tmp]['IS'] = tcounter_reads_indels, tcounter_reads_sf
        vcf_in_row.samples[nid_tmp]['IS'] = ncounter_reads_indels, ncounter_reads_sf
        len_ref = len(ref)
        len_alts = len(alts)
        alt = alts[0]
        if len_alts == 1:
            len_alt = len(alts[0])
        else:
            tgenotype = tcall['GT']
            idx_alt = tgenotype[1]
            alt = alts[idx_alt-1]
            len_alt = len(alt)
        if type == str_cplx or type == str_isrt or type == str_dele or len_ref != len_alt:
            #if (tcounter_reads_indels + tcounter_reads_sf) >= tvd:
            pert_tnoise = float(tcounter_reads_indels + tcounter_reads_sf - tvd) / float(tdp)
            pert_nnoise = float(ncounter_reads_indels + ncounter_reads_sf) / float(ndp)
            # pert_nnoise = float(ncounter_reads_indels + ncounter_reads_sf - nvd) / float(ndp)
            if pert_tnoise > tnoise or pert_nnoise > nnoise:
                vcf_in_row.filter.add('CPX')
                out_forR.append(
                    str(chr) + "\t" + str(pos) + "\t" + "CPX" + "\t" + str(pert_tnoise) + "\t" + str(pert_nnoise))
            else:
                out_forR.append(
                    str(chr) + "\t" + str(pos) + "\t" + "PASS" + "\t" + str(pert_tnoise) + "\t" + str(pert_nnoise))
        vcf_out_fw.write(vcf_in_row)
    tbam.close()
    nbam.close()
    vcf_out_fw.close()
    vcf_in_fr.close()

    # write score results for R
    with open(vcf_out_forR, 'w') as fvcfoutforR:
        # first print column names
        fvcfoutforR.write("CHR" + "\t" + "POS" + "\t" + "TYPE" + "\t" + "T_Noise" + "\t" + "N_Noise" + "\n")
        for tmp_item in out_forR:
            fvcfoutforR.write(tmp_item + "\n")

    # write and run a R script to generate contamination plot in PDF format
    outcpxrdir = os.path.dirname(os.path.realpath(vcf_out_forR))
    writeCPXRScript(outcpxrdir, vcf_out_forR, outr)
    cmd = ["Rscript", outr]
    execute_shell(cmd)


def writeCPXRScript(indir, inf, outf):

    outpdf = outf.replace(".R", ".pdf")
    with open(outf, 'w') as rfile:
        strout = """
        library(ggplot2)
        setwd("%s")
        incpx <- "%s"
        outcpx <- "%s"
        cpx <- read.delim(incpx, header=TRUE)
        pdf(file=outcpx, width=10, height=10)
        ggplot(cpx, aes(x=T_Noise, y=N_Noise, color=TYPE)) + 
        geom_point() +
        labs(x="Tumor Noise",y="Normal Noise",title="Scatter plot for VarDict's Complex/Insertion/Deletion events") +
        theme(plot.title=element_text(size=20,face="bold",vjust=4, hjust=0.5, margin=margin(t=10,b=20)), 
            plot.margin=unit(c(2,2,2,3), "lines"), 
            axis.title=element_text(size=18,face="bold",colour="black"), 
            axis.text=element_text(size=14,face="bold",colour="black"), 
            axis.text.x=element_text(hjust=1, margin=margin(t=20,r=0,b=0,l=0)), 
            axis.text.y=element_text(margin=margin(t=0,r=20,b=0,l=0)), 
            legend.position="bottom", 
            legend.box="horizontal", 
            legend.spacing=unit(1,"cm"), 
            legend.text=element_text(size=14,face="bold",colour="black", margin=margin(r=50,unit="pt")),
            legend.title=element_text(size=14,face="bold",colour="black"),
            legend.key.height=unit(2,"line"), 
            legend.key.width=unit(2,"line"))
        dev.off()
        """
        strout = inspect.cleandoc(strout)
        strout = strout % (indir, inf, outpdf)
        rfile.write(strout)

def execute_shell(cmd):
    try:
        print("Executing %s" % " ".join(cmd), file=sys.stderr)
        subprocess.check_call(" ".join(cmd), shell=True)
    except:
        print("Unexpected Error: %s" % sys.exc_info()[0], file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    print(totaltime)
    sys.exit(0)


