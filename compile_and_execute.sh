#!/bin/bash

#               Author: Matheus S. S. Macedo

#Declared variables
delete_files=false
num_sims=20
num_gens=8000
save_bool=0
save_freq=1000

echo -e "\n\nThis script will compile pgtype_experiment.cpp and run its executable pgtype.exe"

echo -e "\nCompiling"
g++ -I /usr/include/eigen3/ -I ./headers/pgtype_experiment.cpp -o pgtype.exe 

echo -e "Compilation completed!"

echo -e "\nTotal number of simulations on the experiment will be $num_sims
with $num_gens generations each.\n"

./pgtype.exe -nsim $num_sims -gens $num_gens -savepop $save_bool -savefreq $save_freq

if [ "$delete_files" == true ]; then

    echo -e "\n\n Only mean phenotype and mean fit files will be kept after running
            \nDeleting gen* and phen* files."

    rm -rv ./outputs/gen*
    rm -rv ./outputs/phen*
fi

echo -e "\n\nEND OF SCRIPT"
