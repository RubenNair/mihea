#!/bin/bash

loops=30;
flags="";

while getopts pl: flag
do
    case "${flag}" in
        l) loops=${OPTARG};;
        p) flags="${flags} -p";;
    esac
done

echo "Running experiments with ${loops} loops per run!"
### F1-F4 experiments
## GAMBIT_K experiments
echo "[PROGRESS] starting GAMBIT_K experiments for F1-F4!"
./experiments_GAMBIT_2014_recreation.sh -f 1 -f 2 -f 3 -f 4 -o output/test -k -l ${loops} ${flags} & # original experiment settings
./experiments_GAMBIT_2014_recreation.sh -f 1 -f 2 -f 3 -f 4 -o output/test -k -l ${loops} ${flags} -e 1 & # Double population size
./experiments_GAMBIT_2014_recreation.sh -f 1 -f 2 -f 3 -f 4 -o output/test -k -l ${loops} ${flags} -b 3 -e 1 & # Triple population size
./experiments_GAMBIT_2014_recreation.sh -f 1 -f 2 -f 3 -f 4 -o output/test -k -l ${loops} ${flags} -e 2 & # Quadruple population size
wait;
echo "[PROGRESS] GAMBIT_K experiments for F1-F4 done!"

## GAMBIT_R experiments
echo "[PROGRESS] starting GAMBIT_R experiments for F1-F4!"
./experiments_GAMBIT_2014_recreation.sh -f 1 -f 2 -f 3 -f 4 -o output/test -l ${loops} ${flags} & # original experiment settings
./experiments_GAMBIT_2014_recreation.sh -f 1 -f 2 -f 3 -f 4 -o output/test -l ${loops} ${flags} -e 1 & # Double population size
./experiments_GAMBIT_2014_recreation.sh -f 1 -f 2 -f 3 -f 4 -o output/test -l ${loops} ${flags} -b 3 -e 1 & # Triple population size
./experiments_GAMBIT_2014_recreation.sh -f 1 -f 2 -f 3 -f 4 -o output/test -l ${loops} ${flags} -e 2 & # Quadruple population size
wait;
echo "[PROGRESS] GAMBIT_R experiments for F1-F4 done!"

### F5 experiments
## GAMBIT_K experiments
echo "[PROGRESS] starting GAMBIT_K experiments for F5!"
./experiments_GAMBIT_2014_recreation.sh -f 5 -o output/test -k -a 5.0 -a 6.0 -a 7.0 -l ${loops} ${flags} & # Try large values of a to find runs with 30/30
./experiments_GAMBIT_2014_recreation.sh -f 5 -o output/test -k -a 4.0 -a 4.5 -a 5.5 -a 6.5 -l ${loops} ${flags} &
./experiments_GAMBIT_2014_recreation.sh -f 5 -o output/test -k -a 5.0 -a 6.0 -a 7.0 -l ${loops} ${flags} -e 1 & # Same, but double population sizes
./experiments_GAMBIT_2014_recreation.sh -f 5 -o output/test -k -a 4.0 -a 4.5 -a 5.5 -a 6.5 -l ${loops} ${flags} -e 1 &
wait;
echo "[PROGRESS] GAMBIT_K experiments for F5 done!"

## GAMBIT_R experiments
echo "[PROGRESS] starting GAMBIT_R experiments for F5!"
./experiments_GAMBIT_2014_recreation.sh -f 5 -o output/test -a 5.0 -a 6.0 -a 7.0 -l ${loops} ${flags} & # Try large values of a to find runs with 30/30
./experiments_GAMBIT_2014_recreation.sh -f 5 -o output/test -a 4.0 -a 4.5 -a 5.5 -a 6.5 -l ${loops} ${flags} & 
./experiments_GAMBIT_2014_recreation.sh -f 5 -o output/test -a 5.0 -a 6.0 -a 7.0 -l ${loops} ${flags} -e 1 & # Same, but double population sizes
./experiments_GAMBIT_2014_recreation.sh -f 5 -o output/test -a 4.0 -a 4.5 -a 5.5 -a 6.5 -l ${loops} ${flags} -e 1 &
wait;
echo "[PROGRESS] GAMBIT_R experiments for F5 done!"