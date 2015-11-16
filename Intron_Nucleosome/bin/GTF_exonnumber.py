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
python CircosConf.py --input circos.config --output pipe.conf

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid
'''
Chr1    MSU_osa1r7      exon    2903    3268    .       +       .       gene_id "LOC_Os01g01010"; transcript_id "LOC_Os01g01010.1";
Chr1    MSU_osa1r7      exon    3354    3616    .       +       .       gene_id "LOC_Os01g01010"; transcript_id "LOC_Os01g01010.1";
Chr1    MSU_osa1r7      exon    4357    4455    .       +       .       gene_id "LOC_Os01g01010"; transcript_id "LOC_Os01g01010.1";
Chr1    MSU_osa1r7      exon    5457    5560    .       +       .       gene_id "LOC_Os01g01010"; transcript_id "LOC_Os01g01010.1";
Chr1    MSU_osa1r7      exon    7136    7944    .       +       .       gene_id "LOC_Os01g01010"; transcript_id "LOC_Os01g01010.1";
Chr1    MSU_osa1r7      exon    8028    8150    .       +       .       gene_id "LOC_Os01g01010"; transcript_id "LOC_Os01g01010.1";
Chr1    MSU_osa1r7      exon    8232    8320    .       +       .       gene_id "LOC_Os01g01010"; transcript_id "LOC_Os01g01010.1";
Chr1    MSU_osa1r7      exon    8408    8608    .       +       .       gene_id "LOC_Os01g01010"; transcript_id "LOC_Os01g01010.1";
Chr1    MSU_osa1r7      exon    9210    9617    .       +       .       gene_id "LOC_Os01g01010"; transcript_id "LOC_Os01g01010.1";
Chr1    MSU_osa1r7      exon    10104   10187   .       +       .       gene_id "LOC_Os01g01010"; transcript_id "LOC_Os01g01010.1";
Chr1    MSU_osa1r7      exon    10274   10430   .       +       .       gene_id "LOC_Os01g01010"; transcript_id "LOC_Os01g01010.1";
Chr1    MSU_osa1r7      exon    10504   10817   .       +       .       gene_id "LOC_Os01g01010"; transcript_id "LOC_Os01g01010.1";
Chr1    MSU_osa1r7      5UTR    2903    3268    .       +       .       gene_id "LOC_Os01g01010"; transcript_id "LOC_Os01g01010.1";
Chr1    MSU_osa1r7      5UTR    3354    3448    .       +       .       gene_id "LOC_Os01g01010"; transcript_id "LOC_Os01g01010.1";
'''
def readtable(infile):
    data = defaultdict(list)
    exon  = defaultdict(int)
    strand= defaultdict(str)
    r = re.compile(r'gene_id "(.*?)";')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t', line)             
                m = r.search(unit[8])
                gene = m.groups(0)[0] if m else 'NA'
                data[gene].append(unit)
                if unit[2] == 'exon':
                    exon[gene] += 1
                    strand[gene] = unit[6]
    for gene in sorted(data.keys()):
        count = 0
        if strand[gene] == '+':
            for feature in data[gene]:
                if feature[2] == 'exon':
                    count += 1
                    feature[8] += ' exon_number "%s";' %(count)
                newline = '\t'.join(feature)
                print newline
        else:
            for feature in data[gene]:
                if feature[2] == 'exon':
                    count += 1
                    feature[8] += ' exon_number "%s";' %(count)
                newline = '\t'.join(feature)
                print newline

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

