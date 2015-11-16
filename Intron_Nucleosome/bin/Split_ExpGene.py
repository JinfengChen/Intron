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
python Split_ExpGene.py --input ../input/MSU7.gene.exon_number.gtf

Expressed gene: 20478
No Expressed gene: 14901
Low Expressed gene: 4521

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

'''
NB1	NB2	NB3	HEG41	HEG42	HEG43
LOC_Os01g01010	2560	1247	1549	861	1612	949
LOC_Os01g01019	5	18	14	6	9	4
LOC_Os01g01030	182	105	132	82	109	79
'''
def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if line.startswith('LOC'): 
                unit = re.split(r'\t',line)
                count = (float(unit[1]) + float(unit[2]) + float(unit[3]))/3
                data[unit[0]] = count
    return data

'''
Chr1    MSU_osa1r7      exon    2903    3268    .       +       .       gene_id "LOC_Os01g01010"; transcript_id "LOC_Os01g01010.1";
Chr1    MSU_osa1r7      exon    3354    3616    .       +       .       gene_id "LOC_Os01g01010"; transcript_id "LOC_Os01g01010.1";
Chr1    MSU_osa1r7      exon    4357    4455    .       +       .       gene_id "LOC_Os01g01010"; transcript_id "LOC_Os01g01010.1";
Chr1    MSU_osa1r7      exon    5457    5560    .       +       .       gene_id "LOC_Os01g01010"; transcript_id "LOC_Os01g01010.1";
'''
def split(infile, reads):
    nh = defaultdict(int)
    nl = defaultdict(int)
    nn = defaultdict(int)
    prefix = os.path.splitext(infile)[0]
    #print prefix
    r = re.compile(r'gene_id "(.*?)";')
    data = defaultdict(str)
    hfile = open(prefix + '.HighExp.gtf', 'w')
    lfile = open(prefix + '.LowExp.gtf', 'w')
    nfile = open(prefix + '.NoExp.gtf', 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if line.startswith('Chr'): 
                unit = re.split(r'\t',line)
                m = r.search(unit[8])
                gene = m.groups(0)[0] if m else 'NA'
                if reads[gene] > 20:
                    nh[gene] = 1
                    print >> hfile, line
                elif reads[gene] < 2:
                    nn[gene] = 1
                    print >> nfile, line       
                else:
                    nl[gene] = 1
                    print >> lfile, line
    print 'Expressed gene: %s\nNo Expressed gene: %s\nLow Expressed gene: %s' %(str(len(nh.keys())), str(len(nn.keys())), str(len(nl.keys())))
    hfile.close()
    lfile.close()
    nfile.close()

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

    reads = readtable('../input/NB_HEG4.count')
    split(args.input, reads)

if __name__ == '__main__':
    main()

