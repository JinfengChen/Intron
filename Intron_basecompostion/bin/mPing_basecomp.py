#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def usage():
    test="name"
    message='''
python mPing_basecomp.py

    '''
    print message


def base_frequency(dna):
    d = {}
    for base in 'ATCG':
        d[base] = dna.count(base)/float(len(dna))
    return d


def fasta(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        #print record.id
        seq = str(record.seq.upper())
        fastaid[record.id] = seq 
    return fastaid


def basecomposition(fastafile):
    win = 20
    five_prim = defaultdict(list)
    three_prim = defaultdict(list)
    for record in SeqIO.parse(fastafile,"fasta"):
        seq = str(record.seq)
        for i in range(len(seq)/2/win):
            start = i*win + 1
            end = start + win + 1
            seq_win = seq[start:end]
            frq = base_frequency(seq_win)
            frq_AT = frq['A'] + frq['T']
            frq_GC = frq['G'] + frq['C']
            
            #print seq_win, frq_AT
            five_prim[i].append(frq_GC)
        seq_rev = seq[::-1]
        for i in range(len(seq)/2/win):
            start = i*win + 1
            end = start + win + 1
            seq_win = seq_rev[start:end]
            frq = base_frequency(seq_win)
            frq_AT = frq['A'] + frq['T']
            frq_GC = frq['G'] + frq['C']
            three_prim[i].append(frq_GC)          
        #print record.id
    five_prim_avg = defaultdict(float)
    three_prim_avg = defaultdict(float)
    for i in sorted(five_prim.keys(), key=int):
        five_prim_avg[i] = mean(five_prim[i])
    for i in sorted(five_prim.keys(), key=int):
        three_prim_avg[i] = mean(three_prim[i])
    ofile = open('mping.basecomp.distr', 'w')
    base_frq = []
    for i in range(10)[::-1]:
        print i, str(three_prim_avg[i])
        base_frq.append(str(three_prim_avg[i]))
    for i in range(10):
        print i, str(five_prim_avg[i])
        base_frq.append(str(five_prim_avg[i]))
    print >> ofile, '\t'.join(base_frq)
    ofile.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gff')
    parser.add_argument('-b', '--bed')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    
    basecomposition('../input/mping.fa')

if __name__ == '__main__':
    main()


