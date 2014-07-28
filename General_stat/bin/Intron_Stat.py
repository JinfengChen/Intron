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
Python Intron_Stat.py --input ../input/MSU_r7.all.final.full.gff3

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

def median(s):
    i = len(s)
    if not i%2:
        return (s[(i/2)-1]+s[i/2])/2.0
    return s[i/2]

'read full gff, store intron length to dict with gene as key'
def readgff(infile):
    data = defaultdict(lambda : list())
    s = re.compile(r'Parent=(.*)')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t', line)
                m = s.search(unit[8])
                #print '%s\t%s' %(unit[0], unit[4])
                if unit[2] == 'intron' and m:
                    #print 'Intron %s\t%s' %(unit[0], unit[4])
                    gene = m.groups(0)[0]
                    intron_len = int(unit[4])-int(unit[3])+1
                    data[gene].append(intron_len)                    
    return data

'How many intron do every gene have: intron number\tgene number'
def gene_sum(gff):
    t_gene = 39978
    data = defaultdict(int)
    n_gene = len(gff.keys())
    for g in sorted(gff.keys()):
        data[len(gff[g])] = data[len(gff[g])] + 1
    nointron_gene = t_gene - n_gene

    ofile = open ('gene.intron_number.stat', 'w')
    print >> ofile, '%s\t%s' %(0, nointron_gene)
    twenty_intron = 0
    for n_intron in sorted(data.keys(),key=int):
        if n_intron >= 15:
            twenty_intron += n_intron
        else:
            print >> ofile, '%s\t%s' %(n_intron, data[n_intron])
    print >> ofile, '%s\t%s' %('>15', twenty_intron)
    ofile.close()

'intron length: length interval\tintron number'
def intron_len_sum(gff):
    data = defaultdict(int) 
    inter = {
         1:'1-125',
         2:'126-250',
         3:'251-500',
         4:'501-1000',
         5:'1001-2000',
         6:'2001-5000',
         7:'5001-10000',
         8:'>10000'
    }
    all_intron = []
    long_intron = []
    for g in sorted(gff.keys()):
        for l in gff[g]:
            all_intron.append(l)
            if l <= 125:
                data[1] += 1
            elif l <= 250:
                data[2] += 1
            elif l <= 500:
                data[3] += 1
            elif l <= 1000:
                data[4] += 1
            elif l <= 2000:
                data[5] += 1
            elif l <= 5000:
                data[6] += 1
            elif l <= 10000:
                data[7] += 1
            else:
                data[8] += 1
                long_intron.append(len(gff[g]))
    print 'Statistic of intron'
    print '%s\t%s\t%s\t%s' %(min(all_intron), max(all_intron), mean(all_intron), median(all_intron))
    print 'How many introns do these long intron containing gene have?'
    print '%s\t%s\t%s\t%s' %(min(long_intron), max(long_intron), mean(long_intron), median(long_intron))

    ofile = open('intron.length.distr.stat', 'w')
    for r in sorted(data.keys(), key=int):
        print >> ofile, '%s\t%s\t%s' %(r, inter[r], data[r])
 
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

    g_intron = readgff(args.input)
    gene_sum(g_intron)
    intron_len_sum(g_intron)

if __name__ == '__main__':
    main()

