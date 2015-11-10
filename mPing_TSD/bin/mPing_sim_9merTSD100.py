import os
for i in range(1,101):
    #cmd = 'python mPing_sim_9merTSDV2_Random.py --input pictogram/somatic.tsd.matrix --output simulateV2_Random_TSD9mer_somaticMat --number %s' %(i)
    #cmd = 'python mPing_sim_9merTSDV2_Random.py --input pictogram/strain.tsd.matrix --output simulateV2_Random_TSD9mer_strainMat --number %s' %(i)
    #cmd = 'python mPing_sim_9merTSDV2_Random.py --input pictogram/RIL.tsd.matrix --output simulateV2_Random_TSD9mer_rilMat --number %s' %(i)
    cmd = 'python mPing_sim_9merTSDV2_chromatin.py --input pictogram/RIL.tsd.matrix --output simulateV2_Chromain_TSD9mer_rilMat --number %s' %(i)
    os.system(cmd)
