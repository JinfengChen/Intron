python ../bin/ToGFF.py --input 4strains.heter
sort -k1,1 -k3,3n temp.gff > 4strains.heter.gff
python ../bin/ToGFF.py --input 4strains.homo
sort -k1,1 -k3,3n temp.gff > 4strains.homo.gff
python ../bin/ToGFF.py --input 4strains.somatic
sort -k1,1 -k3,3n temp.gff > 4strains.somatic.gff
python ../bin/ToGFF.py --input NB_mping
sort -k1,1 -k3,3n temp.gff > NB_mping.gff
python ../bin/ToGFF.py --input heterozygous.RIL.uniq
sort -k1,1 -k3,3n temp.gff > heterozygous.RIL.uniq.gff
python ../bin/ToGFF.py --input homozygous.RIL.uniq
sort -k1,1 -k3,3n temp.gff > homozygous.RIL.uniq.gff
python ../bin/ToGFF.py --input somatic.RIL.uniq
sort -k1,1 -k3,3n temp.gff > somatic.RIL.uniq.gff

echo "Strains and RILs are mping subjected to selection, somatics are those newly transposition"
cat NB_mping.gff 4strains.homo.gff | sort -k1,1 -k3,3n > Strains.gff
cp heterozygous.RIL.uniq.gff > RIL.gff
cat 4strains.heter.gff 4strains.somatic.gff heterozygous.RIL.uniq.gff somatic.RIL.uniq.gff | sort -k1,1 -k3,3n > Somatic.gff
