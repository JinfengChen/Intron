#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
import math
from scipy import stats
import glob

def usage():
    test="name"
    message='''
python Sim_Sum.py --input simulation

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

'''
input:
1	1-125	35	0.0759219088937
2	126-250	42	0.0911062906725
3	251-500	72	0.156182212581
4	501-1000	115	0.249457700651
5	1001-2000	110	0.238611713666
6	2001-5000	73	0.158351409978
7	5001-10000	9	0.0195227765727
8	>10000	5	0.0108459869848

output:
rank	range	mean	frequency	left_ci95	right_ci95	left_ci95	right_ci95
1	1-125	47.91	0.0928670780684	46.6859662644	49.1340337356	0.0905943113657	0.0951398447712
2	126-250	38.73	0.0751188346511	37.5834267128	39.8765732872	0.0729046487489	0.0773330205533
3	251-500	81.69	0.158278083706	79.7903572904	83.5896427096	0.154901189349	0.161654978064
4	501-1000	141.92	0.275078572028	139.667696129	144.172303871	0.271213751694	0.278943392361
5	1001-2000	118.47	0.229724183309	116.218454179	120.721545821	0.225505273475	0.233943093143
6	2001-5000	67.15	0.130055655723	65.4405196017	68.8594803983	0.127091649949	0.133019661496
7	5001-10000	14.1	0.0273321080793	13.3655655649	14.8344344351	0.0259140544668	0.0287501616917
8	>10000	5.95	0.0115454844356	5.47062751456	6.42937248544	0.0106202923026	0.0124706765686
'''

def readtable(files, outfile):
    ofile = open(outfile, 'w')
    name = defaultdict(lambda: str)
    number = defaultdict(list)
    freq   = defaultdict(list)
    for infile in files:
        filehd = open (infile, 'r')
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                name[unit[0]] = unit[1]
                number[unit[0]].append(float(unit[2]))
                freq[unit[0]].append(float(unit[3]))
        filehd.close()
    for r in sorted(name.keys(), key=int):
        nm, nl, nh = ci95(number[r])
        fm, fl, fh = ci95(freq[r])
        print >> ofile, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(r, name[r], nm, fm, nl, nh, fl, fh)
    ofile.close()


'''
input:
1	35	
2	42
3	72
4		11

output:
rank	mean	left_ci95	right_ci95
1	47.91	46.6859662644	49.1340337356
2	38.73	37.5834267128	39.8765732872
3	81.69	79.7903572904	83.5896427096
4	141.92	139.667696129	144.172303871
'''

def readtable1(files, outfile):
    ofile = open(outfile, 'w')
    name = defaultdict(lambda: str)
    number = defaultdict(list)
    freq   = defaultdict(list)
    for infile in files:
        filehd = open (infile, 'r')
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                number[unit[0]].append(float(unit[1]))
        filehd.close()
    for r in sorted(number.keys(), key=int):
        nm, nl, nh = ci95(number[r])
        print >> ofile, '%s\t%s\t%s\t%s' %(r, nm, nl, nh)
    ofile.close()


'''
intergenic	2474	0.807968647943
five_prime_UTR	95	0.0310254735467
CDS	91	0.0297191378184
intron	281	0.0917700849118
three_prime_UTR	92	0.0300457217505
mRNA	28	0.00914435009798
'''

def readtable2(files, outfile):
    ofile = open(outfile, 'w')
    name = defaultdict(lambda: str)
    number = defaultdict(list)
    freq   = defaultdict(list)
    for infile in files:
        filehd = open (infile, 'r')
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                number[unit[0]].append(float(unit[1]))
                freq[unit[0]].append(float(unit[2]))
        filehd.close()
    for r in sorted(number.keys()):
        nm, nl, nh = ci95(number[r])
        fm, fl, fh = ci95(freq[r])
        print >> ofile, '%s\t%s\t%s\t%s\t%s\t%s\t%s' %(r, nm, fm, nl, nh, fl, fh)
    ofile.close()





'''
sem: standard error of mean
two sem is the 95% confidential interval
http://www.randalolson.com/2012/08/06/statistical-analysis-made-easy-in-python/
'''
def ci95_sem(x):
    mean = np.mean(x)
    sem = stats.sem(x)
    left = mean - 2*sem
    right= mean + 2*sem
    return mean, left, right

def ci95(x):
    mean = np.mean(x)
    sd   = np.std(x) 
    n    = len(x)
    error= 1.96*sd/math.sqrt(n)
    left = mean - error
    right= mean + error
    return mean, left, right

def intron_length(filepath, prefix):
    files = glob.glob(filepath + '/*.intron.length.distr')
    readtable(files, prefix + '.intron.length.distr')

def intron_distance(filepath, prefix):
    files = glob.glob(filepath + '/*.intron.distance.distr')
    readtable(files, prefix + '.intron.distance.distr')

def gene_distance(filepath, prefix):
    files = glob.glob(filepath + '/*.gene.distance.distr')
    readtable(files, prefix + '.gene.distance.distr')    

def intron_position(filepath, prefix):
    files = glob.glob(filepath + '/*.intron_pos.distr')
    readtable1(files, prefix + '.intron_pos.distr') 

def gene_position(filepath, prefix):   
    files = glob.glob(filepath + '/*.position.distr')
    readtable(files, prefix + '.position.distr')

def Five_distance(filepath, prefix):
    files = glob.glob(filepath + '/*.mRNA.5primer.distance.distr')
    readtable(files, prefix + '.mRNA.5primer.distance.distr')

def Three_distance(filepath, prefix):
    files = glob.glob(filepath + '/*.mRNA.3primer.distance.distr')
    readtable(files, prefix + '.mRNA.3primer.distance.distr')

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

    if args.output is None:
        args.output = 'Simulation_TSD9mer_somaticMat'    
 
    intron_length(args.input, args.output)
    intron_distance(args.input, args.output)
    gene_distance(args.input, args.output)
    intron_position(args.input, args.output)
    gene_position(args.input, args.output)
    Five_distance(args.input, args.output)
    Three_distance(args.input, args.output)
 
    x = [2,3,5,6,9]
    ci95(x)
    ci95_sem(x)
if __name__ == '__main__':
    main()

