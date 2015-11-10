#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO
import random

def usage():
    test="name"
    message='''
python mPing_sim_9merTSDV2_AddtionalAnalysis.py --input simulateV2_TSD9mer_rilMat

    '''
    print message


def simulateAA(sim_dir):
    data = defaultdict(list)
    for r in range(1,101):
        sufix = '%04d' %(int(r))
        filename = sim_dir + '/Simulate' + sufix + '.gff'
        intersect_mRNA = os.path.splitext(filename)[0] + '.mRNA.intersect'
        intersect = os.path.splitext(filename)[0] + '.intersect'
        #/rhome/cjinfeng/BigData/00.RD/Intron/Intron_mPing/input/MSU_r7.all.final.mRNA.gff
        #mPing_intergenic.py --input RIL.mRNA.intersect
        #cmd1 = '/opt/bedtools/2.17.0-25-g7b42b3b/bin/bedtools closest -a %s -b /rhome/cjinfeng/BigData/00.RD/Intron/Intron_mPing/input/MSU_r7.all.final.mRNA.gff -d > %s' % (filename, intersect_mRNA)
        #cmd2 = 'python /rhome/cjinfeng/BigData/00.RD/Intron/Intron_mPing/bin/mPing_intergenic.py --input %s' %(intersect_mRNA)
        cmd3 = 'python /rhome/cjinfeng/BigData/00.RD/Intron/Intron_mPing/bin/mPing_position.py --input %s --mrna %s' %(intersect, intersect_mRNA)
        #cmd7 = 'mv %s* %s' %(os.path.splitext(filename)[0], outdir)
        #os.system(cmd1)
        #os.system(cmd2)
        os.system(cmd3)
        #os.system(cmd4)
        #os.system(cmd5)
        #os.system(cmd6)
        #os.system(cmd7)
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()

    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)
   
 
    simulateAA(args.input) 

if __name__ == '__main__':
    main()

