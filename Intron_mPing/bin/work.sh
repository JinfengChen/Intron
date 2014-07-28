bedtools intersect -a ../input/Somatic.gff -b ../input/MSU_r7.all.final.full.utr.gff3 -wao > Somatic.intersect
python mPing_position.py --input Somatic.intersect
python mPing_intron.py --input Somatic.intersect

bedtools intersect -a ../input/Strains.gff -b ../input/MSU_r7.all.final.full.utr.gff3 -wao > Strains.intersect
python mPing_position.py --input Strains.intersect
python mPing_intron.py --input Strains.intersect 


bedtools intersect -a ../input/RIL.gff -b ../input/MSU_r7.all.final.full.utr.gff3 -wao > RIL.intersect
python mPing_position.py --input RIL.intersect
python mPing_intron.py --input RIL.intersect     
