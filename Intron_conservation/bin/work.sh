perl getGene_Os.pl ../input/MSU_r7.all.final.full.gff ../input/MSU_r7.fa --type exon > Os.exon.fa
perl getGene_Ob.pl ../input/Oryza_brachyantha.genome.chr.v1.4.full.gff ../input/Oryza_brachyantha.genome.chr.v1.4.fa --type exon > Ob.exon.fa

blastall -p blastn -i Os.exon.fa -d Ob.exon.fa -e 1e-5 -o Os2Ob.exon.blastm8 -m 8 
echo "one2one ortholog between Os and Ob, 19300"
awk '$3 == 1' ../input/colinear.lst | wc -l


python Intron_cons.py --query ../input/MSU_r7.all.final.full.gff --target ../input/Oryza_brachyantha.genome.chr.v1.4.full.gff > Os_Ob.orth.intron.length
