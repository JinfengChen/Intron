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
perl Intron_cons.py --query ../input/MSU_r7.all.final.full.gff --target ../input/Oryza_brachyantha.genome.chr.v1.4.full.gff

    '''
    print message

'''
#Os_chr10:
Os10g01044.1    Ob10g10010.1    1
Os10g01060.1    Ob10g10030.1    1
Os10g01080.1    Ob10g10050.1    1
'''
def readcol(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'Os'): 
                unit = re.split(r'\t',line)
                if unit[2] == '1':
                    gid = 'LOC_' + unit[0]
                    #print gid, unit[1]
                    data[gid] = unit[1]
    return data

'''
LOC_Os01g01010.1_E1     Ob01g10010.1_E1 94.48   163     9       0       6       168     3       165     6e-67    252
LOC_Os01g01010.1_E2     Ob01g10010.1_E2 94.95   99      5       0       1       99      1       99      1e-38    157
'''
def readblast(infile, collinear):
    data = defaultdict(str)
    p = re.compile(r'(.*)_\w+')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                m1 = p.search(unit[0])
                m2 = p.search(unit[1])
                #print unit[0]
                if m1 and m2:
                    gene1 = m1.groups(0)[0]
                    gene2 = m2.groups(0)[0]
                    #print gene1, gene2, collinear[gene1]
                    if collinear[gene1] == gene2:
                        data[unit[0]] = unit[1]
                        #print unit[0], unit[1]
    return data





'''
Chr1    MSU_osa1r7      mRNA    11218   12435   .       +       .       ID=LOC_Os01g01019.1;Name=expressed protein;
Chr1    MSU_osa1r7      CDS     11798   12060   .       +       .       Parent=LOC_Os01g01019.1
Chr1    MSU_osa1r7      intron  12061   12151   .       +       .       Parent=LOC_Os01g01019.1
Chr1    MSU_osa1r7      CDS     12152   12317   .       +       .       Parent=LOC_Os01g01019.1
'''
def readgff_Os(infile):
    data = defaultdict(list)
    strand = defaultdict(str)
    p = re.compile(r'Parent=(.*)')
    p1 = re.compile(r'ID=(.*?);')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if unit[2] == 'CDS':
                    m = p.search(unit[8])
                    if m:
                        gene = m.groups(0)[0]
                        data[gene].append([unit[3], unit[4]])
                        if not strand.has_key(gene):
                            strand[gene] = unit[6]
    structure1 = defaultdict(list)
    structure2 = defaultdict(list)
    length = defaultdict(int)
    for g in sorted(data.keys()):
        exons = data[g]
        if strand[g] == '+':
            exons.sort(key=lambda x:int(x[0]))
            #print '>',g
            last = exons[0][1]
            last_exon_id = '%s_E%s' %(g, 1)
            for i in range(1,len(exons)):
                #print exons[i][0], exons[i][1]
                intron_id = '%s_I%s' %(g, i)
                intron_len= abs(int(exons[i][0]) - int(last))
                exon_id = '%s_E%s' %(g, str(i+1))
                structure1[intron_id] = [last_exon_id, exon_id]
                structure2[last_exon_id + '_' + exon_id] = intron_id
                length[intron_id] = intron_len
                last = exons[i][1]
                last_exon_id = exon_id
        elif strand[g] == '-':
            exons.sort(key=lambda x:int(x[0]), reverse=True)
            #print '>', g
            last = exons[0][0]
            last_exon_id = '%s_E%s' %(g, 1)
            for i in range(1,len(exons)):
                #print exons[i][0], exons[i][1]
                intron_id = '%s_I%s' %(g, i)
                intron_len= abs(int(exons[i][1]) - int(last))
                exon_id = '%s_E%s' %(g, str(i+1))
                structure1[intron_id] = [last_exon_id, exon_id]
                structure2[last_exon_id + '_' + exon_id] = intron_id
                length[intron_id] = intron_len
                last = exons[i][0]
                last_exon_id = exon_id
    #for intron in sorted(structure1.keys()):
        #print intron, structure1[intron][0], structure1[intron][1]
    #for exon in sorted(structure2.keys()):
        #print exon, structure2[exon]
    return structure1, structure2, length

'''
chr01   Gramene.v1.4    mRNA    88006   89105   .       +       .       ID=Ob01g10080.1;Parent=Ob01g10080;Evidence=Os01t0101200-02:257:0.839869281045752
chr01   Gramene.v1.4    CDS     88006   88041   .       +       0       Parent=Ob01g10080.1;
chr01   Gramene.v1.4    intron  88042   88129   .       +       0       Parent=Ob01g10080.1;
chr01   Gramene.v1.4    CDS     88130   88644   .       +       0       Parent=Ob01g10080.1;
chr01   Gramene.v1.4    intron  88645   88746   .       +       2       Parent=Ob01g10080.1;
chr01   Gramene.v1.4    CDS     88747   88865   .       +       2       Parent=Ob01g10080.1;
chr01   Gramene.v1.4    intron  88866   89007   .       +       1       Parent=Ob01g10080.1;
chr01   Gramene.v1.4    CDS     89008   89105   .       +       1       Parent=Ob01g10080.1;
'''
def readgff_Ob(infile):
    data = defaultdict(list)
    strand = defaultdict(str)
    p = re.compile(r'Parent=(.*);')
    p1 = re.compile(r'ID=(.*?);')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if unit[2] == 'CDS':
                    m = p.search(unit[8])
                    if m:
                        gene = m.groups(0)[0]
                        data[gene].append([unit[3], unit[4]])
                        if not strand.has_key(gene):
                            strand[gene] = unit[6]
    structure1 = defaultdict(list)
    structure2 = defaultdict(list)
    length = defaultdict(int)
    for g in sorted(data.keys()):
        exons = data[g]
        if strand[g] == '+':
            exons.sort(key=lambda x:int(x[0]))
            #print '>',g
            last = exons[0][1]
            last_exon_id = '%s_E%s' %(g, 1)
            for i in range(1,len(exons)):
                #print exons[i][0], exons[i][1]
                intron_id = '%s_I%s' %(g, i)
                intron_len= abs(int(exons[i][0]) - int(last))
                exon_id = '%s_E%s' %(g, str(i+1))
                structure1[intron_id] = [last_exon_id, exon_id]
                structure2[last_exon_id + '_' + exon_id] = intron_id
                length[intron_id] = intron_len
                last = exons[i][1]
                last_exon_id = exon_id
        elif strand[g] == '-':
            exons.sort(key=lambda x:int(x[0]), reverse=True)
            #print '>', g
            last = exons[0][0]
            last_exon_id = '%s_E%s' %(g, 1)
            for i in range(1,len(exons)):
                #print exons[i][0], exons[i][1]
                intron_id = '%s_I%s' %(g, i)
                intron_len= abs(int(exons[i][1]) - int(last))
                exon_id = '%s_E%s' %(g, str(i+1))
                structure1[intron_id] = [last_exon_id, exon_id]
                structure2[last_exon_id + '_' + exon_id] = intron_id
                length[intron_id] = intron_len
                last = exons[i][0]
                last_exon_id = exon_id
    #for intron in sorted(structure1.keys()):
        #print intron, structure1[intron][0], structure1[intron][1]
    #for exon in sorted(structure2.keys()):
        #print exon, structure2[exon]
    return structure1, structure2, length


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--query')
    parser.add_argument('-t', '--target')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.query) > 0
    except:
        usage()
        sys.exit(2)
   
    qry_d1, qry_d2, qry_len = readgff_Os(args.query)
    trg_d1, trg_d2, trg_len = readgff_Ob(args.target)
    collinear = readcol('../input/colinear.lst')
    orthexon = readblast('Os2Ob.exon.blastm8', collinear)    
    
    ofile = open('Os_Ob.orth.intron.length', 'w') 
    for intron in sorted(qry_d1.keys()):
        e1 = qry_d1[intron][0]
        e2 = qry_d1[intron][1]
        #print '>', intron, e1, e2
        if orthexon.has_key(e1) and orthexon.has_key(e2):
            newid = orthexon[e1] + '_' + orthexon[e2]
            if trg_d2.has_key(newid):
                #print orthexon[e1], orthexon[e2], trg_d2[newid]
                intron_ort = trg_d2[newid]
                print >> ofile, '%s\t%s\t%s\t%s' %(intron, qry_len[intron], intron_ort, trg_len[intron_ort])
    ofile.close()

if __name__ == '__main__':
    main()

