perl TSDfrq.pl --input ../input/mping_flank_600.fa > TSD.txt
python mPing_sim.py
python mPing_sim100.py
python Sim_Sum.py --input simulation

cat simulate_intron_length.R | R --slave
cat simulate_intron_distance.R | R --slave

echo "TSD for somatic, RIL and strain"
python GFF2mPing.py --input ../input/Strains.gff > mping_flank_strain.fa
python GFF2mPing.py --input ../input/Somatic.gff > mping_flank_somatic.fa
python GFF2mPing.py --input ../input/RIL.gff > mping_flank_RIL.fa


perl TSD9mer.pl --input ../input/mping_flank_600.fa > 600.tsd.fa 2> 600.tsd.frq
perl TSD9mer.pl --input ../input/mping_flank_somatic.fa > somatic.tsd.fa 2> somatic.tsd.frq
perl TSD9mer.pl --input ../input/mping_flank_strain.fa > strain.tsd.fa 2> strain.tsd.frq
perl TSD9mer.pl --input ../input/mping_flank_RIL.fa > RIL.tsd.fa 2> RIL.tsd.frq

echo "simulate"
python mPing_sim_9merTSD.py --input pictogram/somatic.tsd.matrix --output simulate_TSD9mer_somaticMat

