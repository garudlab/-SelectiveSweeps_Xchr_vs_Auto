#!/bin/bash

num_sims=500000
for NAm in 1110000 15984500; do
    for EUr in 700000; do 
    #401 windows
    hard_in_401=outFiles/Admixture_fixedPopSize_NAm${NAm}_EUr${EUr}_theta_0.01_selection_True_MS_snps_401.txt
    soft_in_401=outFiles/Admixture_fixedPopSize_NAm${NAm}_EUr${EUr}_theta_10_selection_True_MS_snps_401.txt
    outFile_401=outFiles/Admixture_fixedPopSize_NAm${NAm}_EUr${EUr}_BFs_numSim_${num_sims}_401_A.txt
    python compute_Bayes_factors.py $hard_in_401 $soft_in_401 $outFile_401 $num_sims

    #265 windows
    hard_in_265=outFiles/Admixture_fixedPopSize_NAm${NAm}_EUr${EUr}_theta_0.01_selection_True_MS_snps_265.txt
    soft_in_265=outFiles/Admixture_fixedPopSize_NAm${NAm}_EUr${EUr}_theta_10_selection_True_MS_snps_265.txt
    outFile_265=outFiles/Admixture_fixedPopSize_NAm${NAm}_EUr${EUr}_BFs_numSim_${num_sims}_265_A.txt
    python compute_Bayes_factors.py $hard_in_265 $soft_in_265 $outFile_265 $num_sims

    #401 windows 265 downsample
    hard_in=outFiles/Admixture_fixedPopSize_NAm${NAm}_EUr${EUr}_theta_0.01_selection_True_MS_snps.txt
    soft_in=outFiles/Admixture_fixedPopSize_NAm${NAm}_EUr${EUr}_theta_10_selection_True_MS_snps.txt
    outFile=outFiles/Admixture_fixedPopSize_NAm${NAm}_EUr${EUr}_BFs_numSim_${num_sims}_A.txt
    python compute_Bayes_factors.py $hard_in $soft_in $outFile $num_sims
    done
done
