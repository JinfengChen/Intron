#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

def usage():
    test="name"
    message='''

python BlackRegion.py
Concatente 12 chromosomes into one sequence and record the junction regions so we can skip these in simulation.
The concatented genome will help us to avoid biased caused by chromosome length.

    '''
    print message

def fasta(fastafile):
    fasta_seq = defaultdict(str)
    p = re.compile(r'(\d+)')
    for record in SeqIO.parse(fastafile,"fasta"):
        m = p.search(record.id)
        rank = m.groups(0)[0] if m else 0
        fasta_seq[rank] = str(record.seq)
    return fasta_seq


def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = unit[1]
    return data

def blackregion(sequences):
    black = []
    black.append(0)
    count = 0
    genome = ''
    ofile = open('MSU7.Simulation.Chr' ,'w')
    for seq_id in sorted(sequences.keys(), key=int):
        genome += sequences[seq_id]
        length = len(sequences[seq_id])
        temp = black[count] + length
        black.append(temp)
        count += 1
        chrn = 'Chr%s' %(seq_id)
        print >> ofile, '%s\t%s\t%s' %(chrn, black[count-1], black[count])
        start = temp-10
        for i in range(11):
            blackpoint = start + i 
            print blackpoint
    ofile.close()
    newrecord = SeqRecord(Seq(genome),id='Chr1',description="")
    ofile = open ('MSU7.Simulation.Genome.fa','w')
    SeqIO.write(newrecord,ofile,"fasta")
    ofile.close()
    

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()

    ref = '/rhome/cjinfeng/BigData/00.RD/seqlib/MSU_r7.fa'
    fasta_seq = fasta(ref)
    blackregion(fasta_seq)

if __name__ == '__main__':
    main()

