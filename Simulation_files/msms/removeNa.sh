#!/bin/bash

. /u/local/Modules/default/init/modules.sh

module load R

for NAm in 1110000 15984500; do
    for EUr in 700000; do 
    #401 windows
    hard_in_401=outFiles/Admixture_fixedPopSize_NAm${NAm}_EUr${EUr}_theta_0.01_selection_True_MS_snps_401.txt
    soft_in_401=outFiles/Admixture_fixedPopSize_NAm${NAm}_EUr${EUr}_theta_10_selection_True_MS_snps_401.txt

    Rscript RemoveNA.R $hard_in_401
    Rscript RemoveNA.R  $soft_in_401

    #265 windows
    hard_in_265=outFiles/Admixture_fixedPopSize_NAm${NAm}_EUr${EUr}_theta_0.01_selection_True_MS_snps_265.txt
    soft_in_265=outFiles/Admixture_fixedPopSize_NAm${NAm}_EUr${EUr}_theta_10_selection_True_MS_snps_265.txt

    Rscript RemoveNA.R $hard_in_265
    Rscript RemoveNA.R $soft_in_265
  
    #401 windows 265 downsample
    hard_in=outFiles/Admixture_fixedPopSize_NAm${NAm}_EUr${EUr}_theta_0.01_selection_True_MS_snps.txt
    soft_in=outFiles/Admixture_fixedPopSize_NAm${NAm}_EUr${EUr}_theta_10_selection_True_MS_snps.txt

    Rscript RemoveNA.R $hard_in
    Rscript RemoveNA.R $soft_in
    
    done
done