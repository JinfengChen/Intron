#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO

def usage():
    test="name"
    message='''
python mPing_position.py --input Somatic.intersect 

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

def position(pos):
    if len(pos) > 1:
        for i in pos:
            if i == 'mRNA':
                continue
            else:
                return i
    else:
        return pos[0]



'''
Chr1    RelocaTE        mPing   10046640        10046642        .       .       .       Strains=HEG4_2;GT=heterozygous  Chr1    MSU_osa1r7      intergenic
'''
def readtable(infile):
    data = defaultdict(list)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping = unit[0] + '_' + unit[3] + '_' + unit[4]
                pt = unit[11]
                data[mping].append(pt)

    sumx = defaultdict(int)
    tt = 0
    for m in data.keys():
        #print m, data[m]
        p = position(data[m]) 
        sumx[p] += 1      
        tt += 1
    types = ['intergenic', 'five_prime_UTR', 'CDS', 'intron', 'three_prime_UTR', 'mRNA']
    prefix = os.path.splitext(infile)[0]
    ofile = open (prefix + '.position.distr', 'w')
    for p in types:
        print >> ofile, '%s\t%s\t%s' %(p,sumx[p],str(float(sumx[p])/int(tt)))

 

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    readtable(args.input)

if __name__ == '__main__':
    main()

