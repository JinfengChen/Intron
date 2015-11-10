#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
from Bio import SeqIO
sys.path.append('/rhome/cjinfeng/BigData/software/ProgramPython/lib')
from utility import gff_parser, createdir

def usage():
    test="name"
    message='''
python Convert_position_to_one_chromosome.py --chr MSU7.chr.inf --gff GSM655033_Rice_Seedling_DHsites.MSU7.Corrected.gff 

Convert gff position from 12 chromosome to 1 chromosome.
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#Chr1    DHS     DHsites 10537   10645   .       +       .       ID=1;
def convert_gff(infile, pos):
    data = defaultdict(lambda : int())
    ofile = open(re.sub(r'.gff', r'.1chr.gff', infile), 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'Chr'):
                unit = re.split(r'\t', line)
                unit[3] = str(pos[unit[0]] + int(unit[3]))
                unit[4] = str(pos[unit[0]] + int(unit[4]))
                unit[0] = 'Chr1'
                print >> ofile, '\t'.join(unit)
            else:
                print >> ofile, line
    ofile.close()

#Chr01	43270923	16701176	17133774
#Chr10	23207287	8100966	8178267
def read_chr(infile):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'Chr'): 
                unit = re.split(r'\t',line)
                chrs = int(re.sub(r'Chr', r'',unit[0]))
                data[chrs] = int(unit[1])
    last = 0
    pos  = defaultdict(lambda : int())
    for c in sorted(data.keys(), key=int):
        pos['Chr%s' %(c)] = last
        last  += data[c]
        #print 'Chr%s\t%s' %(c, pos[c])
    return pos


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--chr')
    parser.add_argument('-g', '--gff')
    args = parser.parse_args()
    try:
        len(args.chr) > 0
    except:
        usage()
        sys.exit(2)

    pos = read_chr(args.chr)
    convert_gff(args.gff, pos)

if __name__ == '__main__':
    main()

