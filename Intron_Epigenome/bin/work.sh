bedtools intersect -a ../input/MSU_r7.all.final.intron.gff3 -b ../input/H3K9me2.RSEG.bed -wao > MSU_r7.all.final.intron.H3K9me2.intersect
python Intron_Epi.py --input MSU_r7.all.final.intron.H3K9me2.intersect

bedtools intersect -a ../input/MSU_r7.fa.RepeatMasker.out.gff -b ../input/H3K9me2.RSEG.bed -wao > MSU_r7.fa.RepeatMasker.H3K9me2.intersect
python TE_Epi.py --input MSU_r7.fa.RepeatMasker.H3K9me2.intersect

