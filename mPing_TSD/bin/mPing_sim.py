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
python mPing_sim100.py

    '''
    print message

def tsd_in_genome(fastafile):
    data = defaultdict(int)
    ofile1 = open('Potential.mPing.Site.position', 'w')
    ofile2 = open('Potential.mPing.Site.frequency', 'w')
    tsds =['TTA', 'TAA', 'TCA', 'TGA', 'TAG', 'CTA', 'TAC', 'GTA']
    for tsd in tsds:
        p = re.compile(r'%s' % tsd)
        for record in SeqIO.parse(fastafile,"fasta"):
            for m in p.finditer(str(record.seq)):
                print >> ofile1, '%s\t%s\t%s\t%s' %(m.group(), record.id, m.start(), m.end())
                data[tsd] += 1
    for tsd in data.keys():
        print >> ofile2, '%s\t%s' %(tsd, str(data[tsd]))

    ofile1.close()
    ofile2.close()

def simulate(infile, outfile):
    data = defaultdict(list)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if unit[0] == 'TCA' or unit[0] == 'TGA':
                    data['TCA'].append([unit[1], unit[2], unit[3]])
                elif unit[0] == 'TTA' or unit[0] == 'TAA':
                    data['TTA'].append([unit[1], unit[2], unit[3]])
                elif unit[0] == 'TAG' or unit[0] == 'CTA':
                    data['TAG'].append([unit[1], unit[2], unit[3]])
                elif unit[0] == 'TAC' or unit[0] == 'GTA':
                    data['TAC'].append([unit[1], unit[2], unit[3]])

    ofile = open ('temp.gff', 'w')
    tsd = ''
    for i in range(1,3069):
        rn = random.randint(1,10000)
        if rn < 9925:
            tsd = 'TTA'
        elif rn < 9986:
            tsd = 'TCA'
        elif rn < 9993:
            tsd = 'TAG'
        else:
            tsd = 'TAC'
        #print tsd
        index = random.randint(1,len(data[tsd]))
        #index = random.sample(range(len(data[tsd])), 1)[0]
        print >> ofile, '%s\tSimulate\tmPing\t%s\t%s\t.\t.\t.\tID=mPing_%s' %(data[tsd][index][0], data[tsd][index][1], data[tsd][index][2], i)
    ofile.close()
    cmd1 = 'sort -k1,1 -k4,4n %s > %s' %('temp.gff', outfile)
    cmd2 = 'rm temp.gff'
    os.system(cmd1)
    os.system(cmd2)
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()

    if args.output is None:
        args.output = 'Simulate0.gff'

    ref = '/rhome/cjinfeng/BigData/00.RD/seqlib/MSU_r7.fa'
    if not os.path.isfile('Potential.mPing.Site.position'):
        tsd_in_genome(ref) 
    simulate('Potential.mPing.Site.position', args.output)
 
if __name__ == '__main__':
    main()

