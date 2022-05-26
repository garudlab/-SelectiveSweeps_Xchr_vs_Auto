#!/bin/bash

num_sims=52000

 
    #401 windows
    hard_in_401=outFiles_tmp/ConstantNe2657111_theta_0.01_selection_True_Q50_ChrA_snps_401.txt
    soft_in_401=outFiles_tmp/ConstantNe2657111_theta_10_selection_True_Q50_ChrA_snps_401.txt
    outFile_401=outFiles_tmp/ConstantNe2657111_BF_${num_sims}_401_A.txt
    python compute_Bayes_factors.py $hard_in_401 $soft_in_401 $outFile_401 $num_sims

    #265 windows
    hard_in_265=outFiles_tmp/ConstantNe2657111_theta_0.01_selection_True_Q50_ChrA_snps_265.txt
    soft_in_265=outFiles_tmp/ConstantNe2657111_theta_10_selection_True_Q50_ChrA_snps_265.txt
    outFile_265=outFiles_tmp/ConstantNe2657111_BF_${num_sims}_265_A.txt
    python compute_Bayes_factors.py $hard_in_265 $soft_in_265 $outFile_265 $num_sims

    #401 windows 265 downsample
    hard_in=outFiles_tmp/ConstantNe2657111_theta_0.01_selection_True_Q50_ChrA_snps.txt
    soft_in=outFiles_tmp/ConstantNe2657111_theta_10_selection_True_Q50_ChrA_snps.txt
    outFile=outFiles_tmp/ConstantNe2657111_theta_BF_${num_sims}_A.txt
    python compute_Bayes_factors.py $hard_in $soft_in $outFile $num_sims
    

