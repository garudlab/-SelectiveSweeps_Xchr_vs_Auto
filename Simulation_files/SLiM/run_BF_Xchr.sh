#!/bin/bash

num_sims=50000

 
    #401 windows
    hard_in_401=outFiles/ConstantNe2657111_theta_0.01_selection_True_Q50_ChrX_snps_401.txt
    soft_in_401=outFiles/ConstantNe2657111_theta_10_selection_True_Q50_ChrX_snps_401.txt
    outFile_401=outFiles/ConstantNe2657111_BF_${num_sims}_401_X.txt
    python compute_Bayes_factors.py $hard_in_401 $soft_in_401 $outFile_401 $num_sims

    #265 windows
    hard_in_265=outFiles/ConstantNe2657111_theta_0.01_selection_True_Q50_ChrX_snps_265.txt
    soft_in_265=outFiles/ConstantNe2657111_theta_10_selection_True_Q50_ChrX_snps_265.txt
    outFile_265=outFiles/ConstantNe2657111_BF_${num_sims}_265_X.txt
    python compute_Bayes_factors.py $hard_in_265 $soft_in_265 $outFile_265 $num_sims

    #401 windows 265 downsample
    hard_in=outFiles/ConstantNe2657111_theta_0.01_selection_True_Q50_ChrX_snps.txt
    soft_in=outFiles/ConstantNe2657111_theta_10_selection_True_Q50_ChrX_snps.txt
    outFile=outFiles/ConstantNe2657111_theta_BF_${num_sims}_X.txt
    python compute_Bayes_factors.py $hard_in $soft_in $outFile $num_sims
    

