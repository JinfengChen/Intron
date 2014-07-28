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

python Intron_position.py --gff ../../Intron_conservation/input/MSU_r7.all.final.full.gff --bed ../input/Somatic.intersect
    '''
    print message

'''
Chr1    RelocaTE        mPing   10046640        10046642        .       .       .       Strains=HEG4_2;GT=heterozygous  Chr1    MSU_osa1r7      intergenic
'''
def readtable(infile, intron_pos):
    data = defaultdict(list)
    dist = defaultdict(list)
    gdist = defaultdict(list)
    intron_mping = defaultdict(int)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping = unit[0] + '_' + unit[3] + '_' + unit[4]
                pt = unit[11]
                if pt == 'intron':
                    intron = unit[9] + '_' + unit[12] + '_' + unit[13]
                    data[mping] = int(int(unit[13])-int(unit[12]))
                    dist[mping] = min(abs(int(int(unit[13])-int(unit[3]))), int(int(unit[13])-int(unit[4])))
                    if abs(int(unit[13]) - int(unit[12]) + 1) > 2000 and intron_pos.has_key(intron):
                        intron_mping[intron_pos[intron]] += 1
                    elif abs(int(unit[13]) - int(unit[12]) + 1) > 2000:
                        print 'not found'
                elif pt == 'intergenic':
                    gdist[mping] = min(abs(int(int(unit[13])-int(unit[3]))), int(int(unit[13])-int(unit[4])))
    prefix = os.path.splitext(infile)[0]
    ofile = open(prefix + '.intron_pos.distr', 'w')
    for p in sorted(intron_mping, key=int):
        print >> ofile, '%s\t%s' %(p, intron_mping[p])
    return intron_mping
    ofile.close()


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
                        data[gene].append([unit[3], unit[4], unit[0]])
                        if not strand.has_key(gene):
                            strand[gene] = unit[6]
    structure1 = defaultdict(list)
    structure2 = defaultdict(list)
    length = defaultdict(int)
    pos = defaultdict(lambda : int())
    pos_s = defaultdict(lambda : int())
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

                intron_s = int(last) + 1
                intron_e = int(exons[i][0]) - 1
                intron_name = exons[i][2] + '_' + str(intron_s) + '_' + str(intron_e)
                
                last = exons[i][1]
                last_exon_id = exon_id 
                if intron_len > 2000:
                    #print '>',g, '+'
                    #print intron_name
                    pp = i if i < 4 else 4
                    pos[intron_name] = pp
                    pos_s[pp] += 1
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

                intron_e = int(last) - 1
                intron_s = int(exons[i][1]) + 1
                intron_name = exons[i][2] + '_' + str(intron_s) + '_' + str(intron_e)

                last = exons[i][0]
                last_exon_id = exon_id
                if intron_len > 2000:
                    #print '>',g, '-'
                    #print intron_name
                    pp = i if i < 4 else 4
                    pos[intron_name] = pp
                    pos_s[pp] += 1
    #for intron in sorted(structure1.keys()):
        #print intron, structure1[intron][0], structure1[intron][1]
    #for exon in sorted(structure2.keys()):
        #print exon, structure2[exon]
    for i in sorted(pos_s.keys(), key=int):
        print i, pos_s[i]
    return pos


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gff')
    parser.add_argument('-b', '--bed')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.gff) > 0
    except:
        usage()
        sys.exit(2)
   
    pos = readgff_Os(args.gff)
    intron_mping = readtable(args.bed, pos) 

if __name__ == '__main__':
    main()

