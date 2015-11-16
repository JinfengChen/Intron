python GTF_exonnumber.py --input ../input/MSU7.gene.gtf > ../input/MSU7.gene.exon_number.gtf
python Split_ExpGene.py --input ../input/MSU7.gene.exon_number.gtf

samtools view -b -o Chr1.unique.bam ../input/Nucleosome.unique.bam Chr1 &
samtools view -b -o DHS.Chr1.unique.bam ../input/DHS.unique.bam Chr1 &

time python TSSprofile.py > log 2> log2 &
time python TSSprofile_V2.py > log 2> log2 &

echo "Nucleosome TSS"
python TSSprofile1.py ../input/MSU7.gene.exon_number.HighExp.gtf Nucleosome_HighExp_Chr1 > log &
python TSSprofile1.py ../input/MSU7.gene.exon_number.NoExp.gtf Nucleosome_NoExp_Chr1 > log &
python TSSprofile1.py ../input/MSU7.gene.exon_number.gtf Nucleosome_ALLgene_Chr1 > log &


echo "DHS TSS"
python TSSprofile1.py ../input/MSU7.gene.exon_number.HighExp.gtf DHS_HighExp_Chr1 > log &
python TSSprofile1.py ../input/MSU7.gene.exon_number.NoExp.gtf DHS_NoExp_Chr1 > log &
python TSSprofile1.py ../input/MSU7.gene.exon_number.gtf DHS_ALLgene_Chr1 > log &

echo "Nucleosome TSD"
python TSDprofile.py ../input/Somatic.gff Nucleosome_Somatic_Chr1 > log

echo "DHS mPing TSD" 
python TSDprofile.py ../input/Somatic.gff DHS_Somatic_Chr1 > log
python TSDprofile_DHS.py ../input/Somatic.gff DHS_Somatic_Chr1_frag > log

