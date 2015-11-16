import sys
import re
import HTSeq
import numpy
import matplotlib as mpl
mpl.use('pdf')
from matplotlib import pyplot

bamfile = HTSeq.BAM_Reader( "Chr1.unique.bam" )
#bamfile = HTSeq.BAM_Reader( "DHS.Chr1.unique.bam" )
#bamfile = HTSeq.BAM_Reader( "../input/Nucleosome.unique.bam" )
#bamfile = HTSeq.BAM_Reader( "../input/DHS.unique.bam" )
#gtffile = HTSeq.GFF_Reader( "../input/MSU7.gene.exon_number.HighExp.gtf" )
#gtffile = HTSeq.GFF_Reader(sys.argv[1])
#gfffile = '../input/Somatic.gff'
gfffile = sys.argv[1]

halfwinwidth = 1000
fragmentsize = 73
total = 60745783.00/1000000 ## nucleosome
gsize = 372000000
readlen = 36

coverage = HTSeq.GenomicArray( "auto", stranded=False, typecode="i" )
for almnt in bamfile:
   if almnt.aligned:
      #almnt.iv.length = fragmentsize
      if almnt.iv.strand == '+':
          almnt.iv.start = almnt.iv.start + fragmentsize - 10
          almnt.iv.end = almnt.iv.start + 20
          #if not almnt.iv.start < 500:
          #    coverage[ almnt.iv ] += 1
      else:
          almnt.iv.start = almnt.iv.start + readlen - fragmentsize - 10
          almnt.iv.end = almnt.iv.start + 20
      print almnt.iv
      if not almnt.iv.start < 500:
          coverage[ almnt.iv ] += 1

#tsspos = set()
#for feature in gtffile:
#   if feature.type == "exon" and feature.attr["exon_number"] == "1":
#      #print feature.iv.start_d_as_pos.pos
#      if feature.iv.start_d_as_pos.pos > 5000:
#          tsspos.add( feature.iv.start_d_as_pos )

##gff
#Chr1    RelocaTE        mPing   10046640        10046642        .       .       .       Strains=HEG4_2;GT=heterozygous
#Chr1    RelocaTE        mPing   10508960        10508962        .       .       .       Strains=RIL;GT=heterozygous
tsdpos = []
with open (gfffile, 'r') as filehd:
    for line in filehd:
        line = line.rstrip()
        if line.startswith('Chr'): 
            unit = re.split(r'\t',line)
            print unit[0], unit[3], unit[4]
            tsdpos.append([unit[0], unit[3], unit[4]])

profile = numpy.zeros( 2*halfwinwidth, dtype='i' )      
for p in tsdpos:
   window = HTSeq.GenomicInterval( p[0], int(p[1]) - halfwinwidth, int(p[1]) + halfwinwidth, "." )
   wincvg = numpy.fromiter( coverage[window], dtype='i', count=2*halfwinwidth )
   profile += wincvg


pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), profile )  
#pyplot.savefig('Nucleosome_HighExp_Chr1.pdf')
pyplot.savefig('%s.pdf' %(sys.argv[2]))

