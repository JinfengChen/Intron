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
qsub mPing_sim_9merTSD100.sh
python Sim_Sum.py --input simulate_TSD9mer_somaticMat --output Simulate.TSD9mer.SomaticMat
python Sim_Sum.py --input simulate_TSD9mer_strainMat --output Simulate.TSD9mer.StrainMat
python Sim_Sum.py --input simulate_TSD9mer_rilMat --output Simulate.TSD9mer.rilMat

echo "concatente chromosome to do simulation"
python BlackRegion.py > MSU7.Simulation.Blacklist
qsub mPing_sim_9merTSD100.sh

python mPing_sim_9merTSDV2_AddtionalAnalysis.py --input simulateV2_TSD9mer_rilMat
python mPing_sim_9merTSDV2_AddtionalAnalysis.py --input simulateV2_TSD9mer_somaticMat
python mPing_sim_9merTSDV2_AddtionalAnalysis.py --input simulateV2_TSD9mer_strainMat

python Sim_Sum.py --input simulateV2_TSD9mer_somaticMat --output Simulate.TSD9mer.SomaticMat
python Sim_Sum.py --input simulateV2_TSD9mer_strainMat --output Simulate.TSD9mer.StrainMat
python Sim_Sum.py --input simulateV2_TSD9mer_rilMat --output Simulate.TSD9mer.rilMat

echo "Random pick, not using TSD frequency"
python mPing_sim_9merTSD100.py

python mPing_sim_9merTSDV2_AddtionalAnalysis.py --input simulateV2_Random_TSD9mer_somaticMat
python mPing_sim_9merTSDV2_AddtionalAnalysis.py --input simulateV2_Random_TSD9mer_rilMat
python mPing_sim_9merTSDV2_AddtionalAnalysis.py --input simulateV2_Random_TSD9mer_strainMat

python Sim_Sum.py --input simulateV2_Random_TSD9mer_somaticMat --output Simulate.TSD9mer.SomaticMat
python Sim_Sum.py --input simulateV2_Random_TSD9mer_strainMat --output Simulate.TSD9mer.StrainMat
python Sim_Sum.py --input simulateV2_Random_TSD9mer_rilMat --output Simulate.TSD9mer.rilMat

echo "add DHS info"
python Convert_position_to_one_chromosome.py --chr ../input/MSU7.chr.inf --gff ../input/GSM655033_Rice_Seedling_DHsites.MSU7.Corrected.gff

