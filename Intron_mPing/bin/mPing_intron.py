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
python mPing_intron.py --input Somatic.intersect 

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
def intron_len_sum(mping, fname):
    data = defaultdict(int) 
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
    tt = 0
    for m in sorted(mping.keys()):
        l = mping[m]
        if l > 0:
            tt += 1
            if l <= 125:
                data[1] += 1
            elif l <= 250:
                data[2] += 1
            elif l <= 500:
                data[3] += 1
            elif l <= 1000:
                data[4] += 1
            elif l <= 2000:
                data[5] += 1
            elif l <= 5000:
                data[6] += 1
            elif l <= 10000:
                data[7] += 1
            else:
                data[8] += 1

    ofile = open(fname, 'w')
    for r in sorted(data.keys(), key=int):
        print >> ofile, '%s\t%s\t%s\t%s' %(r, inter[r], data[r], float(data[r])/int(tt))


'''
Chr1    RelocaTE        mPing   10046640        10046642        .       .       .       Strains=HEG4_2;GT=heterozygous  Chr1    MSU_osa1r7      intergenic
'''
def readtable(infile):
    data = defaultdict(list)
    dist = defaultdict(list)
    gdist = defaultdict(list)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping = unit[0] + '_' + unit[3] + '_' + unit[4]
                pt = unit[11]
                if pt == 'intron':
                    data[mping] = int(int(unit[13])-int(unit[12]))
                    dist[mping] = min(abs(int(int(unit[13])-int(unit[3]))), int(int(unit[13])-int(unit[4])))
                    #dist[mping] = abs(int(int(unit[13])-int(unit[3])))
                elif pt == 'intergenic':
                    gdist[mping] = min(abs(int(int(unit[13])-int(unit[3]))), int(int(unit[13])-int(unit[4])))
    
    prefix = os.path.splitext(infile)[0]
    intron_len_sum(data, prefix + '.intron.length.distr')
    intron_len_sum(dist, prefix + '.intron.distance.distr')
    intron_len_sum(gdist, prefix + '.gene.distance.distr')
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

