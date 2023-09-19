#!/bin/bash
base=2;
exponent=0;
useKrystofVersion="false";
printGenerationalOutputEachRun="false";
flags="-F2_0_0_-1_0 -N1";
consoleOutput="/dev/stdout";
problems=();
basePopsizes=();
folderBaseName="output";
loops=1;
a_vals=("1.15" "2.75" "3.00")

## Original possible population sizes from 2014 paper (eyeballed from graphs)
ps=("15" "20" "25" "30" "40" "50" "60" "80" "100" "120" "150")

while getopts b:e:kpo:f:l:a: flag
do
    case "${flag}" in
        b) base=${OPTARG};;
        e) exponent=${OPTARG};;
        k) useKrystofVersion="true";;
        p) printGenerationalOutputEachRun="true";;
        o) folderBaseName=${OPTARG};;
        f) problems+=(${OPTARG});;
        l) loops=${OPTARG};;
        a) a_vals+=("${OPTARG}");;
    esac
done

## If a values (relevant for F5) were added to a_vals array, then remove the first three default values and only use specified ones
if [ ${#a_vals[@]} -gt 3 ];
then
    a_vals=("${a_vals[@]:3}");
fi

# echo "${a_vals[@]}"

## Calculate population size multiplier, round to 5 decimals (to avoid some floating point errors)
popsizeMult=$(bc -l <<< "e(${exponent}*l(${base}))");
printf -v popsizeMult %.5f "${popsizeMult}";

# for ((i=1; i<= $loops; i++)); do echo "RUN1 ${i}, PID $!"; done &
# for ((i=1; i<= $loops; i++)); do echo "RUN2 ${i}, PID $!"; done &
# jobs;
# wait;
# (for ((i=1; i<= $loops; i++)); do echo "RUN3 ${i}, PID $!"; done) &
# (for ((i=1; i<= $loops; i++)); do echo "RUN4 ${i}, PID $!"; done) &
# (for ((i=1; i<= $loops; i++)); do echo "RUN5 ${i}, PID $!"; done) &
# (for ((i=1; i<= $loops; i++)); do echo "RUN6 ${i}, PID $!"; done) &
# wait;

if [ "${useKrystofVersion}" = "true" ];
then
    basePopsizes=("20" "30" "100" "200" "150");
    flags="${flags} -k"
    folderBaseName="${folderBaseName}/GAMBIT_K";
else
    basePopsizes=("10" "10" "50" "80" "150")
    folderBaseName="${folderBaseName}/GAMBIT_R";
fi

if [ "${printGenerationalOutputEachRun}" = "false" ];
then
    consoleOutput="/dev/null";
fi

echo -e "[INFO] \tRunning experiment with popsizeMult ${popsizeMult} and Krystof version is ${useKrystofVersion}";
echo -e "[INFO] \tFlags are ${flags}, folderBaseName: ${folderBaseName}";
# build/MixedIntegerGOMEA_O -P1 -n$(echo 15*${popsizeMult} | bc) -L20_20 -V0.0 -O${folderBaseName}/testing $flags & #-k 
# wait;


if [[ "${problems[@]}" =~ 1 ]];
then
    ### COMMANDS FOR F1:
    echo -e "[RUNNING] \tStarting F1 experiments";
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P1 -n$(echo ${basePopsizes[0]}*${popsizeMult} | bc) -L20_20 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F1/0.5/20_20 ${flags} >${consoleOutput}; done & 
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P1 -n$(echo ${basePopsizes[0]}*${popsizeMult} | bc) -L30_10 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F1/0.25/30_10 ${flags} >${consoleOutput}; done & 
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P1 -n$(echo ${basePopsizes[0]}*${popsizeMult} | bc) -L10_30 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F1/0.75/10_30 ${flags} >${consoleOutput}; done & 

    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P1 -n$(echo ${basePopsizes[0]}*${popsizeMult} | bc) -L40_40 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F1/0.5/40_40 ${flags} >${consoleOutput}; done & 
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P1 -n$(echo ${basePopsizes[0]}*${popsizeMult} | bc) -L60_20 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F1/0.25/60_20 ${flags} >${consoleOutput}; done & 
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P1 -n$(echo ${basePopsizes[0]}*${popsizeMult} | bc) -L20_60 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F1/0.75/20_60 ${flags} >${consoleOutput}; done & 

    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P1 -n$(echo ${basePopsizes[0]}*${popsizeMult} | bc) -L60_60 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F1/0.5/60_60 ${flags} >${consoleOutput}; done & 
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P1 -n$(echo ${basePopsizes[0]}*${popsizeMult} | bc) -L90_30 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F1/0.25/90_30 ${flags} >${consoleOutput}; done & 
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P1 -n$(echo ${basePopsizes[0]}*${popsizeMult} | bc) -L30_90 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F1/0.75/30_90 ${flags} >${consoleOutput}; done & 

    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P1 -n$(echo ${basePopsizes[0]}*${popsizeMult} | bc) -L80_80 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F1/0.5/80_80 ${flags} >${consoleOutput}; done & 
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P1 -n$(echo ${basePopsizes[0]}*${popsizeMult} | bc) -L120_40 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F1/0.25/120_40 ${flags} >${consoleOutput}; done & 
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P1 -n$(echo ${basePopsizes[0]}*${popsizeMult} | bc) -L40_120 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F1/0.75/40_120 ${flags} >${consoleOutput}; done & 
    wait;
    echo -e "[RUNNING] \tDone running F1 experiments";
fi

if [[ "${problems[@]}" =~ 2 ]];
then
    ### COMMANDS FOR F2:
    echo -e "[RUNNING] \tStarting F2 experiments";
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P2 -n$(echo ${basePopsizes[1]}*${popsizeMult} | bc) -L20_20 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F2/0.5/20_20 ${flags} >${consoleOutput}; done &
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P2 -n$(echo ${basePopsizes[1]}*${popsizeMult} | bc) -L30_10 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F2/0.25/30_10 ${flags} >${consoleOutput}; done &
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P2 -n$(echo ${basePopsizes[1]}*${popsizeMult} | bc) -L10_30 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F2/0.75/10_30 ${flags} >${consoleOutput}; done &

    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P2 -n$(echo ${basePopsizes[1]}*${popsizeMult} | bc) -L40_40 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F2/0.5/40_40 ${flags} >${consoleOutput}; done &
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P2 -n$(echo ${basePopsizes[1]}*${popsizeMult} | bc) -L60_20 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F2/0.25/60_20 ${flags} >${consoleOutput}; done &
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P2 -n$(echo ${basePopsizes[1]}*${popsizeMult} | bc) -L20_60 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F2/0.75/20_60 ${flags} >${consoleOutput}; done &

    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P2 -n$(echo ${basePopsizes[1]}*${popsizeMult} | bc) -L60_60 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F2/0.5/60_60 ${flags} >${consoleOutput}; done &
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P2 -n$(echo ${basePopsizes[1]}*${popsizeMult} | bc) -L90_30 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F2/0.25/90_30 ${flags} >${consoleOutput}; done &
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P2 -n$(echo ${basePopsizes[1]}*${popsizeMult} | bc) -L30_90 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F2/0.75/30_90 ${flags} >${consoleOutput}; done &

    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P2 -n$(echo ${basePopsizes[1]}*${popsizeMult} | bc) -L80_80 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F2/0.5/80_80 ${flags} >${consoleOutput}; done &
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P2 -n$(echo ${basePopsizes[1]}*${popsizeMult} | bc) -L120_40 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F2/0.25/120_40 ${flags} >${consoleOutput}; done &
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P2 -n$(echo ${basePopsizes[1]}*${popsizeMult} | bc) -L40_120 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F2/0.75/40_120 ${flags} >${consoleOutput}; done &
    wait;
    echo -e "[RUNNING] \tDone running F2 experiments";
fi

if [[ "${problems[@]}" =~ 3 ]];
then
    ### COMMANDS FOR F3:
    echo -e "[RUNNING] \tStarting F3 experiments";
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P3 -n$(echo ${basePopsizes[2]}*${popsizeMult} | bc) -L20_20 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F3/0.5/20_20 ${flags} >${consoleOutput}; done &
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P3 -n$(echo ${basePopsizes[2]}*${popsizeMult} | bc) -L30_10 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F3/0.25/30_10 ${flags} >${consoleOutput}; done &
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P3 -n$(echo ${basePopsizes[2]}*${popsizeMult} | bc) -L10_30 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F3/0.75/10_30 ${flags} >${consoleOutput}; done &

    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P3 -n$(echo ${basePopsizes[2]}*${popsizeMult} | bc) -L40_40 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F3/0.5/40_40 ${flags} >${consoleOutput}; done &
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P3 -n$(echo ${basePopsizes[2]}*${popsizeMult} | bc) -L60_20 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F3/0.25/60_20 ${flags} >${consoleOutput}; done &
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P3 -n$(echo ${basePopsizes[2]}*${popsizeMult} | bc) -L20_60 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F3/0.75/20_60 ${flags} >${consoleOutput}; done &

    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P3 -n$(echo ${basePopsizes[2]}*${popsizeMult} | bc) -L60_60 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F3/0.5/60_60 ${flags} >${consoleOutput}; done &
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P3 -n$(echo ${basePopsizes[2]}*${popsizeMult} | bc) -L90_30 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F3/0.25/90_30 ${flags} >${consoleOutput}; done &
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P3 -n$(echo ${basePopsizes[2]}*${popsizeMult} | bc) -L30_90 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F3/0.75/30_90 ${flags} >${consoleOutput}; done &

    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P3 -n$(echo ${basePopsizes[2]}*${popsizeMult} | bc) -L80_80 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F3/0.5/80_80 ${flags} >${consoleOutput}; done &
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P3 -n$(echo ${basePopsizes[2]}*${popsizeMult} | bc) -L120_40 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F3/0.25/120_40 ${flags} >${consoleOutput}; done &
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P3 -n$(echo ${basePopsizes[2]}*${popsizeMult} | bc) -L40_120 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F3/0.75/40_120 ${flags} >${consoleOutput}; done &
    wait;
    echo -e "[RUNNING] \tDone running F3 experiments";
fi

if [[ "${problems[@]}" =~ 4 ]];
then
    ### COMMANDS FOR F4:
    echo -e "[RUNNING] \tStarting F4 experiments";
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P4 -n$(echo ${basePopsizes[3]}*${popsizeMult} | bc) -L20_20 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F4/0.5/20_20 ${flags} >${consoleOutput}; done &
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P4 -n$(echo ${basePopsizes[3]}*${popsizeMult} | bc) -L30_10 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F4/0.25/30_10 ${flags} >${consoleOutput}; done &
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P4 -n$(echo ${basePopsizes[3]}*${popsizeMult} | bc) -L10_30 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F4/0.75/10_30 ${flags} >${consoleOutput}; done &

    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P4 -n$(echo ${basePopsizes[3]}*${popsizeMult} | bc) -L40_40 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F4/0.5/40_40 ${flags} >${consoleOutput}; done &
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P4 -n$(echo ${basePopsizes[3]}*${popsizeMult} | bc) -L60_20 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F4/0.25/60_20 ${flags} >${consoleOutput}; done &
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P4 -n$(echo ${basePopsizes[3]}*${popsizeMult} | bc) -L20_60 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F4/0.75/20_60 ${flags} >${consoleOutput}; done &

    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P4 -n$(echo ${basePopsizes[3]}*${popsizeMult} | bc) -L60_60 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F4/0.5/60_60 ${flags} >${consoleOutput}; done &
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P4 -n$(echo ${basePopsizes[3]}*${popsizeMult} | bc) -L90_30 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F4/0.25/90_30 ${flags} >${consoleOutput}; done &
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P4 -n$(echo ${basePopsizes[3]}*${popsizeMult} | bc) -L30_90 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F4/0.75/30_90 ${flags} >${consoleOutput}; done &

    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P4 -n$(echo ${basePopsizes[3]}*${popsizeMult} | bc) -L80_80 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F4/0.5/80_80 ${flags} >${consoleOutput}; done &
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P4 -n$(echo ${basePopsizes[3]}*${popsizeMult} | bc) -L120_40 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F4/0.25/120_40 ${flags} >${consoleOutput}; done &
    for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P4 -n$(echo ${basePopsizes[3]}*${popsizeMult} | bc) -L40_120 -V0.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F4/0.75/40_120 ${flags} >${consoleOutput}; done &
    wait;
    echo -e "[RUNNING] \tDone running F4 experiments";
fi

wait;

if [[ "${problems[@]}" =~ 5 ]];
then
    # Loop over a_vals array
    for a_val in "${a_vals[@]}"; 
    do
        echo -e "[RUNNING] \tStarting F5 experiments for a_val = ${a_val}";
        for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P5_${a_val} -n$(echo ${basePopsizes[4]}*${popsizeMult} | bc) -L10_10 -V2.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F5/a_${a_val}/10_10 ${flags} >${consoleOutput}; done &
        for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P5_${a_val} -n$(echo ${basePopsizes[4]}*${popsizeMult} | bc) -L20_20 -V4.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F5/a_${a_val}/20_20 ${flags} >${consoleOutput}; done &
        for ((i=1; i<= $loops; i++)); do timeout 10m build/MixedIntegerGOMEA_O -P5_${a_val} -n$(echo ${basePopsizes[4]}*${popsizeMult} | bc) -L30_30 -V6.0 -O${folderBaseName}/popsizeMult_${popsizeMult}/${loops}_runs/F5/a_${a_val}/30_30 ${flags} >${consoleOutput}; done &
    done
    wait;
    echo -e "[RUNNING] \tDone running all F5 experiments";
fi



############################################################################################################################################
##################################################### ORIGINAL AND ALT. F5 EXPERIMENTS #####################################################
############################################################################################################################################

## ### ORIGINAL COMMANDS FOR F5:
## build/MixedIntegerGOMEA_O -P5_1.15 -n150 -L10_10 -V2.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R/F5_2xpopsize/a_1.15/10_10 -F2_0_0_10_0 -N1 #-k
## build/MixedIntegerGOMEA_O -P5_1.15 -n200 -L20_20 -V4.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R/F5_2xpopsize/a_1.15/20_20 -F2_0_0_20_0 -N1 #-k
## build/MixedIntegerGOMEA_O -P5_1.15 -n250 -L30_30 -V6.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R/F5_2xpopsize/a_1.15/30_30 -F2_0_0_30_0 -N1 #-k

## build/MixedIntegerGOMEA_O -P5_2.75 -n100 -L10_10 -V2.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R/F5_2xpopsize/a_2.75/10_10 -F2_0_0_10_0 -N1 #-k
## build/MixedIntegerGOMEA_O -P5_2.75 -n150 -L20_20 -V4.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R/F5_2xpopsize/a_2.75/20_20 -F2_0_0_20_0 -N1 #-k
## build/MixedIntegerGOMEA_O -P5_2.75 -n200 -L30_30 -V6.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R/F5_2xpopsize/a_2.75/30_30 -F2_0_0_30_0 -N1 #-k

## build/MixedIntegerGOMEA_O -P5_3.0 -n100 -L10_10 -V2.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R/F5_2xpopsize/a_3.0/10_10 -F2_0_0_10_0 -N1 #-k
## build/MixedIntegerGOMEA_O -P5_3.0 -n100 -L20_20 -V4.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R/F5_2xpopsize/a_3.0/20_20 -F2_0_0_20_0 -N1 #-k
## build/MixedIntegerGOMEA_O -P5_3.0 -n120 -L30_30 -V6.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R/F5_2xpopsize/a_3.0/30_30 -F2_0_0_30_0 -N1 #-k


# ### COMMANDS FOR F5 NOT ROTATED!:
# build/MixedIntegerGOMEA_O -P55_1.15 -n220 -L10_10 -V2.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F55_2xpopsize/a_1.15/10_10 -F2_0_0_10_0 -N1 #-k
# build/MixedIntegerGOMEA_O -P55_1.15 -n220 -L20_20 -V4.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F55_2xpopsize/a_1.15/20_20 -F2_0_0_20_0 -N1 #-k
# build/MixedIntegerGOMEA_O -P55_1.15 -n240 -L30_30 -V6.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F55_2xpopsize/a_1.15/30_30 -F2_0_0_30_0 -N1 #-k

# build/MixedIntegerGOMEA_O -P55_2.75 -n200 -L10_10 -V2.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F55_2xpopsize/a_2.75/10_10 -F2_0_0_10_0 -N1 #-k
# build/MixedIntegerGOMEA_O -P55_2.75 -n220 -L20_20 -V4.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F55_2xpopsize/a_2.75/20_20 -F2_0_0_20_0 -N1 #-k
# build/MixedIntegerGOMEA_O -P55_2.75 -n220 -L30_30 -V6.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F55_2xpopsize/a_2.75/30_30 -F2_0_0_30_0 -N1 #-k

# build/MixedIntegerGOMEA_O -P55_3.0 -n200 -L10_10 -V2.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F55_2xpopsize/a_3.0/10_10 -F2_0_0_10_0 -N1 #-k
# build/MixedIntegerGOMEA_O -P55_3.0 -n200 -L20_20 -V4.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F55_2xpopsize/a_3.0/20_20 -F2_0_0_20_0 -N1 #-k
# build/MixedIntegerGOMEA_O -P55_3.0 -n220 -L30_30 -V6.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F55_2xpopsize/a_3.0/30_30 -F2_0_0_30_0 -N1 #-k

# ### COMMANDS FOR F5 FULLY ROTATED, WRONG EXPONENT!:
# build/MixedIntegerGOMEA_O -P51_1.15 -n220 -L10_10 -V2.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F51_2xpopsize/a_1.15/10_10 -F2_0_0_10_0 -N1 #-k
# build/MixedIntegerGOMEA_O -P51_1.15 -n220 -L20_20 -V4.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F51_2xpopsize/a_1.15/20_20 -F2_0_0_20_0 -N1 #-k
# build/MixedIntegerGOMEA_O -P51_1.15 -n240 -L30_30 -V6.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F51_2xpopsize/a_1.15/30_30 -F2_0_0_30_0 -N1 #-k

# build/MixedIntegerGOMEA_O -P51_2.75 -n200 -L10_10 -V2.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F51_2xpopsize/a_2.75/10_10 -F2_0_0_10_0 -N1 #-k
# build/MixedIntegerGOMEA_O -P51_2.75 -n220 -L20_20 -V4.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F51_2xpopsize/a_2.75/20_20 -F2_0_0_20_0 -N1 #-k
# build/MixedIntegerGOMEA_O -P51_2.75 -n220 -L30_30 -V6.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F51_2xpopsize/a_2.75/30_30 -F2_0_0_30_0 -N1 #-k

# build/MixedIntegerGOMEA_O -P51_3.0 -n200 -L10_10 -V2.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F51_2xpopsize/a_3.0/10_10 -F2_0_0_10_0 -N1 #-k
# build/MixedIntegerGOMEA_O -P51_3.0 -n200 -L20_20 -V4.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F51_2xpopsize/a_3.0/20_20 -F2_0_0_20_0 -N1 #-k
# build/MixedIntegerGOMEA_O -P51_3.0 -n220 -L30_30 -V6.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F51_2xpopsize/a_3.0/30_30 -F2_0_0_30_0 -N1 #-k

# ### COMMANDS FOR F5 NOT ROTATED, WRONG EXPONENT!:
# build/MixedIntegerGOMEA_O -P551_1.15 -n220 -L10_10 -V2.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F551_2xpopsize/a_1.15/10_10 -F2_0_0_10_0 -N1 #-k
# build/MixedIntegerGOMEA_O -P551_1.15 -n220 -L20_20 -V4.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F551_2xpopsize/a_1.15/20_20 -F2_0_0_20_0 -N1 #-k
# build/MixedIntegerGOMEA_O -P551_1.15 -n240 -L30_30 -V6.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F551_2xpopsize/a_1.15/30_30 -F2_0_0_30_0 -N1 #-k

# build/MixedIntegerGOMEA_O -P551_2.75 -n200 -L10_10 -V2.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F551_2xpopsize/a_2.75/10_10 -F2_0_0_10_0 -N1 #-k
# build/MixedIntegerGOMEA_O -P551_2.75 -n220 -L20_20 -V4.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F551_2xpopsize/a_2.75/20_20 -F2_0_0_20_0 -N1 #-k
# build/MixedIntegerGOMEA_O -P551_2.75 -n220 -L30_30 -V6.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F551_2xpopsize/a_2.75/30_30 -F2_0_0_30_0 -N1 #-k

# build/MixedIntegerGOMEA_O -P551_3.0 -n200 -L10_10 -V2.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F551_2xpopsize/a_3.0/10_10 -F2_0_0_10_0 -N1 #-k
# build/MixedIntegerGOMEA_O -P551_3.0 -n200 -L20_20 -V4.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F551_2xpopsize/a_3.0/20_20 -F2_0_0_20_0 -N1 #-k
# build/MixedIntegerGOMEA_O -P551_3.0 -n220 -L30_30 -V6.0 -Ooutput/after_update_meeting_31_07_GAMBIT_R_init_out_of_optimum_range/F551_2xpopsize/a_3.0/30_30 -F2_0_0_30_0 -N1 #-k

# ### TEST FOR PETER-> COMMANDS FOR F5, CENTERS OF ELLIPSOIDS ALL AT 0:
# build/MixedIntegerGOMEA_O -P50_1.15 -n220 -L10_10 -V2.0 -Ooutput/after_update_meeting_31_07_GAMBIT_K/F50_doublepopsize/a_1.15/10_10 -F2_0_0_10_0 -N1 -k
# build/MixedIntegerGOMEA_O -P50_1.15 -n220 -L20_20 -V4.0 -Ooutput/after_update_meeting_31_07_GAMBIT_K/F50_doublepopsize/a_1.15/20_20 -F2_0_0_20_0 -N1 -k
# build/MixedIntegerGOMEA_O -P50_1.15 -n240 -L30_30 -V6.0 -Ooutput/after_update_meeting_31_07_GAMBIT_K/F50_doublepopsize/a_1.15/30_30 -F2_0_0_30_0 -N1 -k

# build/MixedIntegerGOMEA_O -P50_2.75 -n200 -L10_10 -V2.0 -Ooutput/after_update_meeting_31_07_GAMBIT_K/F50_doublepopsize/a_2.75/10_10 -F2_0_0_10_0 -N1 -k
# build/MixedIntegerGOMEA_O -P50_2.75 -n220 -L20_20 -V4.0 -Ooutput/after_update_meeting_31_07_GAMBIT_K/F50_doublepopsize/a_2.75/20_20 -F2_0_0_20_0 -N1 -k
# build/MixedIntegerGOMEA_O -P50_2.75 -n220 -L30_30 -V6.0 -Ooutput/after_update_meeting_31_07_GAMBIT_K/F50_doublepopsize/a_2.75/30_30 -F2_0_0_30_0 -N1 -k

# build/MixedIntegerGOMEA_O -P50_3.0 -n200 -L10_10 -V2.0 -Ooutput/after_update_meeting_31_07_GAMBIT_K/F50_doublepopsize/a_3.0/10_10 -F2_0_0_10_0 -N1 -k
# build/MixedIntegerGOMEA_O -P50_3.0 -n200 -L20_20 -V4.0 -Ooutput/after_update_meeting_31_07_GAMBIT_K/F50_doublepopsize/a_3.0/20_20 -F2_0_0_20_0 -N1 -k
# build/MixedIntegerGOMEA_O -P50_3.0 -n220 -L30_30 -V6.0 -Ooutput/after_update_meeting_31_07_GAMBIT_K/F50_doublepopsize/a_3.0/30_30 -F2_0_0_30_0 -N1 -k
