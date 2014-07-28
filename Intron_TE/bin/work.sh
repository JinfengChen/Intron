bedtools intersect -a ../input/MSU_r7.all.final.intron.gff3 -b ../input/MSU_r7.fa.RepeatMasker.out.gff -wao > MSU_r7.all.final.intron.TE.intersect
python Intron_TE.py --input MSU_r7.all.final.intron.TE.intersect

