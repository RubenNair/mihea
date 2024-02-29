#!/bin/bash

loops=30;
flags="";
repeats=1;

while getopts pl:r: flag
do
    case "${flag}" in
        l) loops=${OPTARG};;
        p) flags="${flags} -p";;
        r) repeats=${OPTARG};;
    esac
done
for ((i=1; i<= $repeats; i++))
do
    outputdir="output/test4/repeat_${i}"
    echo "Running experiments with ${loops} loops per run!"
    ### F1-F4 experiments
    ## GAMBIT_K experiments
    echo "[PROGRESS] starting GAMBIT_K experiments for F1-F4!"
    ./experiments_GAMBIT_2014_recreation.sh -f 1 -f 2 -f 3 -f 4 -o ${outputdir} -k -l ${loops} ${flags} -b 2 -e 0 &
    ./experiments_GAMBIT_2014_recreation.sh -f 1 -f 2 -f 3 -f 4 -o ${outputdir} -k -l ${loops} ${flags} -b 2 -e 1 &
    ./experiments_GAMBIT_2014_recreation.sh -f 1 -f 2 -f 3 -f 4 -o ${outputdir} -k -l ${loops} ${flags} -b 2 -e 2 &
    ./experiments_GAMBIT_2014_recreation.sh -f 1 -f 2 -f 3 -f 4 -o ${outputdir} -k -l ${loops} ${flags} -b 2 -e 3 &
    wait;
    echo "[PROGRESS] GAMBIT_K experiments for F1-F4 done!"

    ## GAMBIT_R experiments
    echo "[PROGRESS] starting GAMBIT_R experiments for F1-F4!"
    ./experiments_GAMBIT_2014_recreation.sh -f 1 -f 2 -f 3 -f 4 -o ${outputdir} -l ${loops} ${flags} -b 2 -e 0 &
    ./experiments_GAMBIT_2014_recreation.sh -f 1 -f 2 -f 3 -f 4 -o ${outputdir} -l ${loops} ${flags} -b 2 -e 1 &
    ./experiments_GAMBIT_2014_recreation.sh -f 1 -f 2 -f 3 -f 4 -o ${outputdir} -l ${loops} ${flags} -b 2 -e 2 &
    ./experiments_GAMBIT_2014_recreation.sh -f 1 -f 2 -f 3 -f 4 -o ${outputdir} -l ${loops} ${flags} -b 2 -e 3 &
    wait;
    echo "[PROGRESS] GAMBIT_R experiments for F1-F4 done!"

    ### F5 experiments
    ## GAMBIT_K experiments
    echo "[PROGRESS] starting GAMBIT_K experiments for F5!"
    ./experiments_GAMBIT_2014_recreation.sh -f 5 -o ${outputdir} -k -a 4.625 -a 4.6875 -a 4.75 -a 4.8125  -a 4.875 -a 4.9375 -a 5.0 -a 5.0625 -a 5.125 -a 5.1875 -a 5.25 -a 5.5 -l ${loops} ${flags} -b 2 -e 0 &
    wait;
    ./experiments_GAMBIT_2014_recreation.sh -f 5 -o ${outputdir} -k -a 4.625 -a 4.6875 -a 4.75 -a 4.8125  -a 4.875 -a 4.9375 -a 5.0 -a 5.0625 -a 5.125 -a 5.1875 -a 5.25 -a 5.5 -l ${loops} ${flags} -b 2 -e 1 &
    wait;
    ./experiments_GAMBIT_2014_recreation.sh -f 5 -o ${outputdir} -k -a 4.625 -a 4.6875 -a 4.75 -a 4.8125  -a 4.875 -a 4.9375 -a 5.0 -a 5.0625 -a 5.125 -a 5.1875 -a 5.25 -a 5.5 -l ${loops} ${flags} -b 2 -e 2 &
    wait;
    ./experiments_GAMBIT_2014_recreation.sh -f 5 -o ${outputdir} -k -a 4.625 -a 4.6875 -a 4.75 -a 4.8125  -a 4.875 -a 4.9375 -a 5.0 -a 5.0625 -a 5.125 -a 5.1875 -a 5.25 -a 5.5 -l ${loops} ${flags} -b 2 -e 3 &
    wait;
    ./experiments_GAMBIT_2014_recreation.sh -f 5 -o ${outputdir} -k -a 4.625 -a 4.6875 -a 4.75 -a 4.8125  -a 4.875 -a 4.9375 -a 5.0 -a 5.0625 -a 5.125 -a 5.1875 -a 5.25 -a 5.5 -l ${loops} ${flags} -b 2 -e 4 &
    wait;
    echo "[PROGRESS] GAMBIT_K experiments for F5 done!"

    ## GAMBIT_R experiments
    echo "[PROGRESS] starting GAMBIT_R experiments for F5!"
    ./experiments_GAMBIT_2014_recreation.sh -f 5 -o ${outputdir} -a 4.20 -a 4.225 -a 4.25 -a 4.28125 -a 4.3125 -a 4.34375 -a 4.375 -a 4.4375 -a 4.5 -l ${loops} ${flags} -b 2 -e 0 &
    wait;
    ./experiments_GAMBIT_2014_recreation.sh -f 5 -o ${outputdir} -a 4.20 -a 4.225 -a 4.25 -a 4.28125 -a 4.3125 -a 4.34375 -a 4.375 -a 4.4375 -a 4.5 -l ${loops} ${flags} -b 2 -e 1 &
    wait;
    ./experiments_GAMBIT_2014_recreation.sh -f 5 -o ${outputdir} -a 4.20 -a 4.225 -a 4.25 -a 4.28125 -a 4.3125 -a 4.34375 -a 4.375 -a 4.4375 -a 4.5 -l ${loops} ${flags} -b 2 -e 2 &
    wait;
    ./experiments_GAMBIT_2014_recreation.sh -f 5 -o ${outputdir} -a 4.20 -a 4.225 -a 4.25 -a 4.28125 -a 4.3125 -a 4.34375 -a 4.375 -a 4.4375 -a 4.5 -l ${loops} ${flags} -b 2 -e 3 &
    wait;
    ./experiments_GAMBIT_2014_recreation.sh -f 5 -o ${outputdir} -a 4.20 -a 4.225 -a 4.25 -a 4.28125 -a 4.3125 -a 4.34375 -a 4.375 -a 4.4375 -a 4.5 -l ${loops} ${flags} -b 2 -e 4 &
    wait;
    echo "[PROGRESS] GAMBIT_R experiments for F5 done!"
done

