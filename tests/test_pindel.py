'''
@Description : This tool helps to test PINDEL
@Created :  03/23/2017
@Updated : 03/23/2017
@author : Ronak H Shah

'''

import unittest
import filecmp
import os
from subprocess import Popen
import shlex
def main():
    this_dir, this_filename = os.path.split(__file__)
    new_dir = os.path.dirname(this_dir)
    inputFileVcf = os.path.join(new_dir, "data", "sample_input", "PoolTumor2-T_bc52_PINDEL_0.2.5a7.vcf")
    outFileVcf = os.path.join(new_dir, "PoolTumor2-T_bc52_PINDEL_0.2.5a7_STDfilter.vcf")
    outFileTxt = os.path.join(new_dir, "PoolTumor2-T_bc52_PINDEL_0.2.5a7_STDfilter.txt")
    cmpFileVcf = os.path.join(new_dir, "data", "sample_output", "PoolTumor2-T_bc52_PINDEL_0.2.5a7_STDfilter.vcf")
    cmpFileTxt = os.path.join(new_dir, "data", "sample_output", "PoolTumor2-T_bc52_PINDEL_0.2.5a7_STDfilter.txt")
    scriptFile = os.path.join(new_dir, "filter_pindel.py")
    cmd = "python " + scriptFile + " -tsn PoolTumor2-T " + "-ivcf " + inputFileVcf
    args = shlex.split(cmd)
    if(os.path.isfile(outFileTxt) or (os.path.isfile(outFileVcf))):
        os.remove(outFileTxt)
        os.remove(outFileVcf)
    else:
        proc = Popen(args)
        proc.wait()
        retcode = proc.returncode
        if(retcode >= 0):
            code = 1
        else:
            assert 0
        assert filecmp.cmp(outFileTxt, cmpFileTxt)
        assert filecmp.cmp(outFileVcf, cmpFileVcf)

if __name__ == '__main__':
    unittest.main()
