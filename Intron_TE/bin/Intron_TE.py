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
python Intron_Epi.py --input MSU_r7.all.final.intron.H3K9me2.intersect

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

def position(pos):
    if len(pos) > 1:
        for i in pos:
            if i == 'mRNA':
                continue
            else:
                return i
    else:
        return pos[0]

''
def intron_len_class(intron):
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
    l1 = intron
    cls = 0
    if l1 > 0: # no use, if intron have length larger than 0
        if l1 <= 125:
            cls = 1
        elif l1 <= 250:
            cls = 2
        elif l1 <= 500:
            cls = 3
        elif l1 <= 1000:
            cls = 4
        elif l1 <= 2000:
            cls = 5
        elif l1 <= 5000:
            cls = 6
        elif l1 <= 10000:
            cls = 7
        else:
            cls = 8
    return cls


'''
Chr1    MSU_osa1r7      intron  5561    7135    .       +       .       Parent=LOC_Os01g01010.1 Chr1    RepeatMasker    Transposon      6949    6996    232     -       .       ID=TE000006;Target=MDM2 213 260;Class=DNA/MuDR;PercDiv=27.1
'''
def readtable(infile):
    introns = defaultdict(int)
    total = defaultdict(lambda : defaultdict(lambda : int()))
    data  = defaultdict(lambda : defaultdict(lambda : int()))
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
    p = re.compile(r'Class=(.*?);')
    CACTA = re.compile(r'Spm', re.IGNORECASE)
    MUDR = re.compile(r'MuDR', re.IGNORECASE)
    MITE = re.compile(r'Stowaway|Tourist|MITE', re.IGNORECASE)
    GYPSY = re.compile(r'Gypsy', re.IGNORECASE)
    COPIA = re.compile(r'Copia', re.IGNORECASE)
    LTR = re.compile(r'^LTR$', re.IGNORECASE)
    DNA = re.compile(r'^DNA$', re.IGNORECASE)
    LINE = re.compile(r'^LINE', re.IGNORECASE)
    SINE = re.compile(r'^SINE', re.IGNORECASE)
    HAT = re.compile(r'hAT', re.IGNORECASE)
    HEL = re.compile(r'helitron', re.IGNORECASE) 
    SAT = re.compile(r'Simple_repeat', re.IGNORECASE)
    header = defaultdict(str)
    family = ['Gypsy', 'Copia', 'LINE', 'SINE', 'CACTA', 'hAT', 'MuDR', 'MITE', 'Unknown']
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                intron = unit[0] + '_' + unit[3] + '_' + unit[4]
                length = int(int(unit[4])-int(unit[3]))
                cls = intron_len_class(length)
                overlap = int(unit[18])
                te_type = 'Unknown'
                match = p.search(unit[17])
                if match:
                    te_class = match.groups(0)[0]
                    if MITE.search(te_class):
                        te_type = 'MITE'
                    elif GYPSY.search(te_class):
                        te_type = 'Gypsy'
                    elif COPIA.search(te_class):
                        te_type = 'Copia'
                    elif MUDR.search(te_class):
                        te_type = 'MuDR'
                    elif CACTA.search(te_class):
                        te_type = 'CACTA'
                        #print line
                        #print te_type, cls ,overlap
                    elif LINE.search(te_class):
                        te_type = 'LINE'
                    elif SINE.search(te_class):
                        te_type = 'SINE'
                    elif HAT.search(te_class):
                        te_type = 'hAT'
                    elif HEL.search(te_class):
                        te_type = 'Helitron'
                    elif SAT.search(te_class):
                        te_type = 'Simple_repeat'
                    elif DNA.search(te_class):
                        te_type = 'DNA_other'
                    elif LTR.search(te_class):
                        te_type = 'LTR_other'
                    else:
                        te_type = 'Unknown'
                if not introns.has_key(intron):
                    total[cls]['total'] += length
                introns[intron] = 1
                data[cls][te_type] += overlap
                header[te_type] =1
    h = "Rank\tRange\tTotal\t" + "\t".join(sorted(header.keys()))
    h1= "Rank\tRange\tTotal\t" + "\t".join(family)
    prefix = os.path.splitext(infile)[0]
    ofile = open (prefix + '.length.TE.class', 'w')
    print >> ofile, h1
    for r in sorted(data.keys(), key=int):
        array = [str(r), inter[r], str(total[r]['total'])]
        for t in family:
            #print r, t, data[r][t]
            array.append(str(data[r][t]))
        line = '\t'.join(array)
        print >> ofile, line
    ofile.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
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

