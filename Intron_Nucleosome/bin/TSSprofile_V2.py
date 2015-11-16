import HTSeq
import numpy
import matplotlib as mpl
mpl.use('pdf')
from matplotlib import pyplot

#bamfile_nucls = HTSeq.BAM_Reader( "Chr1.unique.bam" )
bamfile_nucls = HTSeq.BAM_Reader( "../input/Nucleosome.unique.bam" )
bamfile_DHS   = HTSeq.BAM_Reader("../input/DHS.unique.bam")
#bamfile = HTSeq.BAM_Reader( "../input/Nucleosome.unique.bam" )
gtffile = HTSeq.GFF_Reader( "../input/MSU7.gene.exon_number.HighExp.gtf" )
halfwinwidth = 1000
fragmentsize = 200
total = 60745783.00/1000000 ## nucleosome
gsize = 372000000

#Nucleosome
coverage_nucls = HTSeq.GenomicArray( "auto", stranded=False, typecode="i" )
for almnt in bamfile_nucls:
   if almnt.aligned:
      #almnt.iv.length = fragmentsize
      #print almnt.iv
      if not almnt.iv.start < 500:
          coverage_nucls[ almnt.iv ] += 1

#DHS
coverage_DHS = HTSeq.GenomicArray( "auto", stranded=False, typecode="i" )
for almnt in bamfile_DHS:
   if almnt.aligned:
      #almnt.iv.length = fragmentsize
      #print almnt.iv
      if not almnt.iv.start < 500:
          coverage_DHS[ almnt.iv ] += 1



tsspos = set()
for feature in gtffile:
   if feature.type == "exon" and feature.attr["exon_number"] == "1":
      #print feature.iv.start_d_as_pos.pos
      if feature.iv.start_d_as_pos.pos > 5000:
          tsspos.add( feature.iv.start_d_as_pos )

profile_nucls = numpy.zeros( 2*halfwinwidth, dtype='i' )
profile_DHS   = numpy.zeros( 2*halfwinwidth, dtype='i' ) 
for p in tsspos:
   window = HTSeq.GenomicInterval( p.chrom, p.pos - halfwinwidth, p.pos + halfwinwidth, "." )
   wincvg_nucls = numpy.fromiter( coverage_nucls[window], dtype='i', count=2*halfwinwidth )
   wincvg_DHS   = numpy.fromiter( coverage_DHS[window], dtype='i', count=2*halfwinwidth )
   if p.strand == "+":
      profile_nucls += wincvg_nucls
      profile_DHS += wincvg_DHS
   else:
      profile_nucls += wincvg_nucls[::-1]
      profile_DHS   += wincvg_DHS[::-1]

for i in range(len(profile_nucls)):
    profile_nucls[i] = profile_nucls[i]/(gsize*total)
    profile_DHS[i] = profile_DHS[i]/(gsize*total)

pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), profile_nucls ) 
pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), profile_DHS ) 
pyplot.savefig('example02.pdf')


