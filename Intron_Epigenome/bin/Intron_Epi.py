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
python Intron_Epi.py --input MSU_r7.all.final.intron.H3K9me2.intersect

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

'intron length: length interval\tintron number'
def intron_len_sum(intron, fname):
    total = defaultdict(int) 
    epi   = defaultdict(int)
    inter = {
         1:'1-125',
         2:'126-250',
         3:'251-500',
         4:'501-1000',
         5:'1001-2000',
         6:'2001-5000',
         7:'5001-10000',
         8:'>10000'
    }
    for i in sorted(intron.keys()):
        l1 = intron[i][0]
        l2 = intron[i][1]
        if l1 > 0: # no use, if intron have length larger than 0
            if l1 <= 125:
                total[1] += l1
                epi[1] += l2
            elif l1 <= 250:
                total[2] += l1
                epi[2] += l2
            elif l1 <= 500:
                total[3] += l1
                epi[3] += l2
            elif l1 <= 1000:
                total[4] += l1
                epi[4] += l2
            elif l1 <= 2000:
                total[5] += l1
                epi[5] += l2
            elif l1 <= 5000:
                total[6] += l1
                epi[6] += l2
            elif l1 <= 10000:
                total[7] += l1
                epi[7] += l2
            else:
                total[8] += l1
                epi[8] += l2

    ofile = open(fname, 'w')
    for r in sorted(total.keys(), key=int):
        print >> ofile, '%s\t%s\t%s' %(r, inter[r], float(epi[r])/float(total[r]))


'''
Chr1    MSU_osa1r7      intron  69963   70026   .       +       .       Parent=LOC_Os01g01140.1 Chr1    66500   72450   ENRICHED        65.8235 16.9571 +       64
'''
def readtable(infile):
    data = defaultdict(list)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                intron = unit[0] + '_' + unit[3] + '_' + unit[4]
                length = int(int(unit[4])-int(unit[3]))
                overlap = int(unit[16])
                data[intron] = [length, overlap]
    prefix = os.path.splitext(infile)[0]
    intron_len_sum(data, prefix + '.distr')

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
    
    readtable(args.input)

if __name__ == '__main__':
    main()

