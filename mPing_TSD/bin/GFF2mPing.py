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
python GFF2mPing.py --input ../input/Strains.gff > mping_flank_strain.fa

    '''
    print message

def fasta(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        #print record.id
        seq = str(record.seq.upper())
        fastaid[record.id] = seq 
    return fastaid



'''
IN:
Chr1    RelocaTE        mPing   10029365        10029367        .       .       .       Strains=RIL;GT=homozygous

OUT:
>mping.Chr10_18900305_18900307 mping.fwd mping A119_2 Chr10:18900305..18900307 FLANK1:1..600 TSD1:598..600 TE:601..1030 TSD2:1031..1033 FLANK2:1031..1630
CACGTTAACGCCACGTGGAAAGAAGACCGAGTCAATACTGTCACGTAGGCGTCACGTCAGCGAAACCACTCTTCAAAATCGTCCAGGGAGTCAAATTGCACCGGTTTTAAGAGTTTGGGAGTCGAGATATCCGGTTTTATGGTTTTTATGGTTTAGGGATACGAATTAGATTTCGATCACTTTTAAGGGTCATGAAGTGAACTTATTCCTTGCAGGGGCT

'''
def readtable(infile, ref):
    mping = 'GGCCAGTCACAATGGGGGTTTCACTGGTGTGTCATGCACATTTAATAGGGGTAAGACTGAATAAAAAATGATTATTTGCATGAAATGGGGATGAGAGAGAAGGAAAGAGTTTCATCCTGGTGAAACTCGTCAGCGTCGTTTCCAAGTCCTCGGTAACAGAGTGAAACCCCCGTTGAGGCCGATTCGTTTCATTCACCGGATCTCTTGCGTCCGCCTCCGCCGTGCGACCTCCGCATTCTCCCGCGCCGCGCCGGATTTTGGGTACAAATGATCCCAGCAACTTGTATCAATTAAATGCTTTGCTTAGTCTTGGAAACGTCAAAGTGAAACCCCTCCACTGTGGGGATTGTTTCATAAAAGATTTCATTTGAGAGAAGATGGTATAATATTTTGGGTAGCCGTGCAATGACACTAGCCATTGTGACTGGCC'
    data = defaultdict(str)
    p = re.compile(r'Strains=(\w+);')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t', line) 
                chrseq = ref[unit[0]]
                start1 = int(unit[4]) - 600
                end1   = int(unit[4])
                start2 = int(unit[3]) - 1
                end2   = int(unit[3]) + 600 -1 
                flank1 = chrseq[start1:end1] 
                flank2 = chrseq[start2:end2]
                m = p.search(unit[8])
                strain = m.groups(0)[0] if m else 'NA'
                header = 'mping.%s_%s_%s mping.fwd mping %s %s:%s..%s FLANK1:1..600 TSD1:598..600 TE:601..1030 TSD2:1031..1033 FLANK2:1031..1630' %(unit[0], unit[3], unit[4], strain, unit[0], unit[3], unit[4])
                sequence = flank1 + mping + flank2
                print '>%s\n%s' %(header, sequence)
    return data


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

    ref = fasta('/rhome/cjinfeng/BigData/00.RD/seqlib/MSU_r7.fa')
    readtable(args.input, ref)

if __name__ == '__main__':
    main()

