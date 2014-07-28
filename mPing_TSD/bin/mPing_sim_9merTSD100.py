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
python mPing_sim_9merTSD.py --input pictogram/somatic.tsd.matrix --output simulate_TSD9mer_somaticMat

    '''
    print message

def fasta(fastafile):
    fastaid = defaultdict(str)
    fastalen = defaultdict(int)
    for record in SeqIO.parse(fastafile,"fasta"):
        #print record.id
        seq = str(record.seq.upper())
        fastaid[record.id] = seq
        fastalen[record.id] = len(seq) 
    return fastaid, fastalen

def complement(s): 
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    letters = list(s) 
    letters = [basecomplement[base] for base in letters] 
    return ''.join(letters)

def revcom(s):
    return complement(s[::-1])


'''
read frquency matirx
>A	C	G	T
0.110	0.335	0.139	0.416
0.247	0.213	0.138	0.402
0.406	0.370	0.144	0.081
0.000	0.000	0.001	0.999
0.478	0.012	0.010	0.500
0.999	0.001	0.000	0.000
0.092	0.149	0.357	0.403
0.397	0.147	0.205	0.252
0.386	0.150	0.346	0.118
'''
def readmatrix(infile):
    data = defaultdict(lambda : defaultdict(lambda: float))
    count = 0
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'0'): 
                unit = re.split(r'\t',line)
                data['A'][count] = float(unit[0])
                data['C'][count] = float(unit[1])
                data['G'][count] = float(unit[2])
                data['T'][count] = float(unit[3])
                count += 1
    return data

def validdna(dna):
    flag = 1
    bases = ['A', 'T', 'C', 'G'] 
    for b in dna:
        if b not in bases:
            flag = 0
    return flag

def insertion_in_genome(fastaseq, fastalen, matrix):
    data = []
    mping = 0
    while (mping == 0):
        chrn = random.randint(1,12)
        chri = 'Chr' + str(chrn)
        seq = fastaseq[chri]
        rpos = random.randint(1,fastalen[chri]-10)
        s    = rpos-1
        e    = s + 9 
        sitep = seq[s:e]
        if not validdna(sitep):
            continue
        siten = revcom(sitep) 
        probp = 1.00
        probn = 1.00
        for i in range(0,9):
            basep = sitep[i]
            basen = siten[i]
            probp = probp*matrix[basep][i]
            probn = probn*matrix[basen][i]
        prob = max(probp, probn)
        site = sitep if probp == prob else siten
        rn = random.random()
        if rn <= prob:
            mping = 1
            start = s + 3
            end   = s + 5
            data = [chri, str(start), str(end), site]
            print '%s\t%s\t%s' %(chri, str(start), str(end))
    return data


def simulate(fastaseq, fastalen, matrix, outdir):
    data = defaultdict(list)
    cmd0 = 'mkdir %s' %(outdir)
    os.system(cmd0)
    for r in range(1,3):
        sufix = '%04d' %(r)
        filename = 'Simulate' + sufix + '.gff'
        ofile = open ('temp.gff', 'w')
        tsd = ''
        for i in range(1,3069):
            mping = insertion_in_genome(fastaseq, fastalen, matrix) 
            print >> ofile, '%s\tSimulate\tmPing\t%s\t%s\t.\t.\t.\tID=mPing_%s' %(mping[0], mping[1], mping[2], i)
        ofile.close()
        cmd1 = 'sort -k1,1 -k4,4n %s > %s' %('temp.gff', filename)
        cmd2 = 'rm temp.gff'
        intersect = os.path.splitext(filename)[0] + '.intersect'
        cmd3 = 'bedtools intersect -a %s -b ../../Intron_mPing/input/MSU_r7.all.final.full.utr.gff3 -wao > %s' % (filename, intersect)
        cmd4 = 'python ../../Intron_mPing/bin/mPing_intron.py --input %s' %(intersect)
        cmd5 = 'python ../../Intron_mPing/bin/mPing_position.py --input %s' %(intersect)
        cmd6 = 'python ../../Intron_position/bin/Intron_position.py --gff ../../Intron_conservation/input/MSU_r7.all.final.full.gff --bed %s' %(intersect)
        cmd7 = 'mv %s* ./simulation' %(os.path.splitext(filename)[0])
        os.system(cmd1)
        os.system(cmd2)
        os.system(cmd3)
        os.system(cmd4)
        os.system(cmd5)
        os.system(cmd6)
        os.system(cmd7)
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output')
    parser.add_argument('-i', '--input')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()

    if args.output is None:
        args.output = 'simulate_TSD9mer'

    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)
   
 
    ref = '/rhome/cjinfeng/BigData/00.RD/seqlib/MSU_r7.fa'
    fastaseq, fastalen = fasta(ref)
    matrix = readmatrix(args.input)
    simulate(fastaseq, fastalen, matrix, args.output)     
    

if __name__ == '__main__':
    main()

