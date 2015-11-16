import HTSeq
import numpy
import matplotlib as mpl
mpl.use('pdf')
from matplotlib import pyplot

bamfile = HTSeq.BAM_Reader( "Chr1.unique.bam" )
#bamfile = HTSeq.BAM_Reader( "../input/Nucleosome.unique.bam" )
gtffile = HTSeq.GFF_Reader( "../input/MSU7.gene.exon_number.HighExp.gtf" )

halfwinwidth = 1000
fragmentsize = 73
total = 60745783.00/1000000 ## nucleosome
gsize = 20478.00*2000.00
readlen = 36

coverage = HTSeq.GenomicArray( "auto", stranded=False, typecode="i" )
for almnt in bamfile:
   if almnt.aligned:
      print almnt.iv
      #almnt.iv.length = fragmentsize 
      if almnt.iv.strand == '+':
          almnt.iv.start = almnt.iv.start + fragmentsize
          almnt.iv.end = almnt.iv.start + 1
          if almnt.iv.start > 500 and almnt.iv.strand == '+':
              coverage[ almnt.iv ] += 1
      else:
          almnt.iv.start = almnt.iv.start + readlen - fragmentsize
          almnt.iv.end = almnt.iv.start + 1
      print almnt.iv

#      if almnt.iv.start > 500 and almnt.iv.strand == '+':
#          coverage[ almnt.iv ] += 1

tsspos = set()
for feature in gtffile:
   if feature.type == "exon" and feature.attr["exon_number"] == "1":
      #print feature.iv.start_d_as_pos.pos
      if feature.iv.start_d_as_pos.pos > 5000:
          tsspos.add( feature.iv.start_d_as_pos )

profile = numpy.zeros( 2*halfwinwidth, dtype=float )      
for p in tsspos:
   window = HTSeq.GenomicInterval( p.chrom, p.pos - halfwinwidth, p.pos + halfwinwidth, "." )
   wincvg = numpy.fromiter( coverage[window], dtype='i', count=2*halfwinwidth )
   if p.strand == "+":
      profile += wincvg
   else:
      profile += wincvg[::-1]



ofile = open('Nucleosome.profile', 'w')
for i in range(len(profile)):
    profile[i] = float(profile[i])/(gsize*total)
    print >> ofile, profile[i]
ofile.close()

pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), profile )
#pyplot.ylim((0,0.00002)) 
pyplot.savefig('example01_chr1_1bp.pdf')

