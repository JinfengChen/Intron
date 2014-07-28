#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO
from subprocess import PIPE, Popen

def usage():
    test="name"
    message='''
python ToGFF.py --input 4strains.heter
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

def cmdline(command):
    process = Popen(
        args=command,
        stdout=PIPE,
        shell=True
    )
    return process.communicate()[0]


'''
A123_2	Chr8:19850940..19850942	heterozygous
Chr1    SV      Deletion        100145  103713  .       .       .       Size=3569;
'''
def readtable(infile):
    outfile = infile + '.gff'
    ofile = open ('temp.gff', 'w')
    s = re.compile(r'(\w+):(\d+)\.\.(\d+)')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                m = s.search(unit[1])
                if m:
                    print >> ofile, '%s\tRelocaTE\tmPing\t%s\t%s\t.\t.\t.\tStrains=%s;GT=%s' %(m.groups(0)[0], m.groups(0)[1], m.groups(0)[2], unit[0], unit[2]) 
    ofile.close
    cmd = 'sort -k1,1 -k3,3n temp.gff > ' + outfile
    print cmd

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

