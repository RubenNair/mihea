#!/bin/bash
trap "echo '[KILL] SIGINT received, killing process.'; exit" INT

printOutput="false";
flags="-F2_0_0_-1_0 -N5 -b0.99";
maxTime="7200";
consoleOutput="/dev/stdout";
folderBaseName="out/experiments/sample_size_scaling_4_nodes_test_small";
loops=15;
dataInputFile="randomnetwork_ew_4nodes";
dataInputRun="0";
pInd=("10400200" "10400300" "10401000" "10402000" "10404000" "10408000" "10416000"); # Problem Indices
sampleSizes=(200 300 1000 2000 4000 8000 16000);
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
        # Run the 3 methods once for each problem
        # Method 1 init 1
        timeout 140m build/MixedIntegerGOMEA_O -P${pInd[$i]} -n100 -I${dataInputFile}_${sampleSizes[$i]}samples -O${folderBaseName}/problem${pInd[$i]}/Method_1_init_1/run${run} ${flags} -r${run} >${consoleOutput} &
        # Method 3 init 1
        timeout 140m build/MixedIntegerGOMEA_O -P${pInd[$i]} -n100 -I${dataInputFile}_${sampleSizes[$i]}samples -O${folderBaseName}/problem${pInd[$i]}/Method_3_init_1/run${run} ${flags} -r${run} -Y >${consoleOutput} &
        # Method 3 init 3
        timeout 140m build/MixedIntegerGOMEA_O -P${pInd[$i]} -n100 -I${dataInputFile}_${sampleSizes[$i]}samples -O${folderBaseName}/problem${pInd[$i]}/Method_3_init_3/run${run} ${flags} -r${run} -Y -g >${consoleOutput} &
        # CBN-GOMEA
        # Wait for processes to finish every 5 runs (also check it is not first run to make sure 5 runs happen before waiting)
        if (( $run > 0 )) && (( $run % 5 == 0 )); then
            wait;
        fi
    done
    wait;
done
