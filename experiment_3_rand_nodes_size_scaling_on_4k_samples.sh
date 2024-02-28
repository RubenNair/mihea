#!/bin/bash
trap "echo '[KILL] SIGINT received, killing process.'; exit" INT

printOutput="false";
flags="-F2_0_0_-1_0 -N1 -b0.99";
maxTime="18000";
consoleOutput="/dev/stdout";
folderBaseName="out/experiments/rand_nodes_scaling_4k_samples";
loops=30;
dataInputFile="randomnetwork_rand";
dataInputRun="0";
pInd=("30304000" "30404000" "30504000"); # Problem Indices
sampleSizes=(4000 4000 4000);
numNodes=(3 4 5);
popSizes=(45 55 64);
maxGens=20;


while getopts po:l:i:m:T:P:  flag
do
    case "${flag}" in
        p) printOutput="true";;
        o) folderBaseName=${OPTARG};;
        l) loops=${OPTARG};;
        i) dataInputFile=${OPTARG};;
        r) dataInputRun=${OPTARG};;
        m) maxGens=${OPTARG};;
        T) maxTime=${OPTARG};;
        P) pInd=${OPTARG};;
    esac
done

# Add max gens, data input file and run to flags
flags="${flags} -T${maxTime}";

if [ "${printOutput}" = "false" ];
then
    consoleOutput="/dev/null";
fi

echo -e "[INFO] \tStarting experiment...";
echo -e "[INFO] \tFlags are ${flags}, folderBaseName: ${folderBaseName}";


for i in "${!pInd[@]}";
do
    # Repeat experiment $loops times
    for ((run=0; run<$loops; run++)); 
    do
        # Wait for processes to finish every 5 runs (also check it is not first run to make sure 5 runs happen before waiting)
        if (( $run > 0 )) && (( $run % 10 == 0 )); then
            wait;
        fi
        # Run the 3 methods once for each problem
        # Method 1 init 1
        timeout 320m build/MixedIntegerGOMEA_O -P${pInd[$i]} -n${popSizes[$i]} -I${dataInputFile}_${numNodes[$i]}nodes_${sampleSizes[$i]}samples -O${folderBaseName}/problem${pInd[$i]}/Method_1_init_1/run${run} ${flags} -r${run} >${consoleOutput} &
        # Method 3 init 1
        timeout 320m build/MixedIntegerGOMEA_O -P${pInd[$i]} -n${popSizes[$i]} -I${dataInputFile}_${numNodes[$i]}nodes_${sampleSizes[$i]}samples -O${folderBaseName}/problem${pInd[$i]}/Method_3_init_1/run${run} ${flags} -r${run} -Y >${consoleOutput} &
        # Method 3 init 3
        timeout 320m build/MixedIntegerGOMEA_O -P${pInd[$i]} -n${popSizes[$i]} -I${dataInputFile}_${numNodes[$i]}nodes_${sampleSizes[$i]}samples -O${folderBaseName}/problem${pInd[$i]}/Method_3_init_3/run${run} ${flags} -r${run} -Y -g >${consoleOutput} &
        # CBN-GOMEA
        # All relevant commands for CBN-GOMEA can be created by running the following command in the Implementations folder of other (own modified version of examine-code by Damy) project:
        #  python run_script.py -d 6 1 1 2 1 3 True 10400000
        # To execute these commands automatically as well, run it without dry-run flag (-d)
    done
    wait;
done
