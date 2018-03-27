'''
@Description : This tool helps to test VarDict 
@Created :  03/23/2017
@Updated : 03/23/2017
@author : Ronak H Shah

'''

import filecmp
import os
from subprocess import Popen
import shlex
import nose
import logging

def setup_module(): 
    this_dir, this_filename = os.path.split(__file__)
    new_dir = os.path.dirname(this_dir)
    inputFileVcf = os.path.join(new_dir, "data", "sample_input", "PoolTumor2-T_bc52_VarDict_1.4.6.vcf")
    outFileVcf = os.path.join(new_dir, "PoolTumor2-T_bc52_VarDict_1.4.6_STDfilter.vcf")
    outFileTxt = os.path.join(new_dir, "PoolTumor2-T_bc52_VarDict_1.4.6_STDfilter.txt")
    cmpFileVcf = os.path.join(new_dir, "data", "sample_output", "PoolTumor2-T_bc52_VarDict_1.4.6_STDfilter.vcf")
    cmpFileTxt = os.path.join(new_dir, "data", "sample_output", "PoolTumor2-T_bc52_VarDict_1.4.6_STDfilter.txt")
    hspFileVcf = os.path.join(new_dir, "data", "hotspot-list-union-v1-v2.vcf")
    scriptFile = os.path.join(new_dir, "filter_vardict.py")
    cmd = "python " + scriptFile + " -v -tsn PoolTumor2-T " + "-ivcf " + inputFileVcf + " -hvcf " + hspFileVcf
    args = shlex.split(cmd)
    if(os.path.isfile(outFileTxt) or (os.path.isfile(outFileVcf))):
        os.remove(outFileTxt)
        os.remove(outFileVcf)
        
    try:
        proc = Popen(args)
        proc.wait()
        retcode = proc.returncode
        if(retcode >= 0):
            pass
    except:
        e = sys.exc_info()[0]
        logging.info("Running of python command: %s \n has failed. The exception produced is %s Thus we will exit",cmd,e)
        sys.exit(1)
             
def teardown_module():
    this_dir, this_filename = os.path.split(__file__)
    new_dir = os.path.dirname(this_dir)
    outFileVcf = os.path.join(new_dir, "PoolTumor2-T_bc52_VarDict_1.4.6_STDfilter.vcf")
    outFileTxt = os.path.join(new_dir, "PoolTumor2-T_bc52_VarDict_1.4.6_STDfilter.txt")
    cmpFileVcf = os.path.join(new_dir, "data", "sample_output", "PoolTumor2-T_bc52_VarDict_1.4.6_STDfilter.vcf")
    cmpFileTxt = os.path.join(new_dir, "data", "sample_output", "PoolTumor2-T_bc52_VarDict_1.4.6_STDfilter.txt")
    if(os.path.isfile(outFileTxt) or (os.path.isfile(outFileVcf))):
        os.remove(outFileTxt)
        os.remove(outFileVcf)

def test_text_fileSimilarity():
    this_dir, this_filename = os.path.split(__file__)
    new_dir = os.path.dirname(this_dir)
    outFileTxt = os.path.join(new_dir, "PoolTumor2-T_bc52_VarDict_1.4.6_STDfilter.txt")
    cmpFileTxt = os.path.join(new_dir, "data", "sample_output", "PoolTumor2-T_bc52_VarDict_1.4.6_STDfilter.txt")
    nose.tools.ok_(filecmp.cmp(outFileTxt, cmpFileTxt), msg="The current result text file and the original result text file for vardict are not the same") 

def test_vcf_fileSimilarity():
    this_dir, this_filename = os.path.split(__file__)
    new_dir = os.path.dirname(this_dir)
    outFileVcf = os.path.join(new_dir, "PoolTumor2-T_bc52_VarDict_1.4.6_STDfilter.vcf")
    cmpFileVcf = os.path.join(new_dir, "data", "sample_output", "PoolTumor2-T_bc52_VarDict_1.4.6_STDfilter.vcf")
    nose.tools.ok_(filecmp.cmp(outFileVcf, cmpFileVcf), msg="The current result vcf file and the original result vcf file for vardict are not the same")

if __name__ == '__main__':
    nose.main()
