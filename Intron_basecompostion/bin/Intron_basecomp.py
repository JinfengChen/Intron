#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def usage():
    test="name"
    message='''
python Intron_basecomp.py --gff ../input/MSU_r7.all.final.full.gff

    '''
    print message


def base_frequency(dna):
    d = {}
    for base in 'ATCG':
        d[base] = dna.count(base)/float(len(dna))
    return d


def fasta(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        #print record.id
        seq = str(record.seq.upper())
        fastaid[record.id] = seq 
    return fastaid

'''
Chr1    MSU_osa1r7      mRNA    11218   12435   .       +       .       ID=LOC_Os01g01019.1;Name=expressed protein;
Chr1    MSU_osa1r7      CDS     11798   12060   .       +       .       Parent=LOC_Os01g01019.1
Chr1    MSU_osa1r7      intron  12061   12151   .       +       .       Parent=LOC_Os01g01019.1
Chr1    MSU_osa1r7      CDS     12152   12317   .       +       .       Parent=LOC_Os01g01019.1
'''
def readgff_Os(infile, ref):
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
    ofile = open('intron.fa', 'w')
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
                temp_seq = ref[exons[i][2]]
                intron_seq = temp_seq[intron_s-1:intron_e]
                #print exons[i][2], intron_seq,intron_s, intron_e
                print >> ofile, '>%s\n%s' %(intron_name, intron_seq) 

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
                temp_seq = ref[exons[i][2]]
                intron_seq = temp_seq[intron_s-1:intron_e]
                my_seq = Seq(intron_seq, IUPAC.unambiguous_dna)
                intron_seq_rc = str(my_seq.reverse_complement())
                print >> ofile, '>%s\n%s' %(intron_name, intron_seq_rc)

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
    ofile.close()
    #for i in sorted(pos_s.keys(), key=int):
    #    print i, pos_s[i]
    return pos

def basecomposition(fastafile, outfile):
    win = 20
    five_prim = defaultdict(list)
    three_prim = defaultdict(list)
    for record in SeqIO.parse(fastafile,"fasta"):
        seq = str(record.seq)
        for i in range(len(seq)/2/win):
            start = i*win + 1
            end = start + win + 1
            seq_win = seq[start:end]
            frq = base_frequency(seq_win)
            frq_AT = frq['A'] + frq['T']
            frq_GC = frq['G'] + frq['C']
            
            #print seq_win, frq_AT
            five_prim[i].append(frq_GC)
        seq_rev = seq[::-1]
        for i in range(len(seq)/2/win):
            start = i*win + 1
            end = start + win + 1
            seq_win = seq_rev[start:end]
            frq = base_frequency(seq_win)
            frq_AT = frq['A'] + frq['T']
            frq_GC = frq['G'] + frq['C']
            three_prim[i].append(frq_GC)          
        #print record.id
    five_prim_avg = defaultdict(float)
    three_prim_avg = defaultdict(float)
    for i in sorted(five_prim.keys(), key=int):
        five_prim_avg[i] = mean(five_prim[i])
    for i in sorted(five_prim.keys(), key=int):
        three_prim_avg[i] = mean(three_prim[i])
    ofile = open(outfile, 'w')
    base_frq = []
    for i in range(100)[::-1]:
        print i, str(three_prim_avg[i])
        base_frq.append(str(three_prim_avg[i]))
    for i in range(100):
        print i, str(five_prim_avg[i])
        base_frq.append(str(five_prim_avg[i]))
    print >> ofile, '\t'.join(base_frq)
    ofile.close()


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
  
    ref = fasta('/rhome/cjinfeng/BigData/00.RD/seqlib/MSU_r7.fa')
    if not os.path.isfile('intron.fa'): 
        pos = readgff_Os(args.gff, ref)
    basecomposition('intron.fa', 'intron.basecomp.distr')

if __name__ == '__main__':
    main()

