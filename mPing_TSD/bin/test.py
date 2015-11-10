#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO
import random


def blackchroms(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                print unit[1], unit[2], unit[0]
                data[unit[1] + '_' + unit[2]] = unit[0] 
    return data

def convert_position(start, end, chroms):
    for p in sorted(chroms.keys()):
        p1 = re.split(r'_', p)   
        print p, p1[0], p1[1], start, end
        if start >= int(p1[0]) and end <= int(p1[1]):
            start1 = start - int(p1[0])
            end1   = end - int(p1[0])
            print chroms[p], start1, end1
            return chroms[p], start1, end1



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output')
    parser.add_argument('-n', '--number')
    parser.add_argument('-i', '--input')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()

    if args.output is None:
        args.output = 'simulate_TSD9mer'

    if args.number is None:
        args.number = '0'

   
 
    chroms = blackchroms('MSU7.Simulation.Chr')
    convert_position(1100,1103, chroms) 

if __name__ == '__main__':
    main()


