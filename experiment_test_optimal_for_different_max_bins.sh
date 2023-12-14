#!/bin/bash
trap "echo '[KILL] SIGINT received, killing process.'; exit" INT

printOutput="false";
flags="-F2_0_0_-1_0 -N1 -b0.99";
consoleOutput="/dev/stdout";
folderBaseName="output/boundary_test_optimal_for_different_max_bins";
loops=10;
dataInputFile="network7";
dataInputRun="0";
maxGens=20;


while getopts po:l:i:m:  flag
do
    case "${flag}" in
        p) printOutput="true";;
        o) folderBaseName=${OPTARG};;
        l) loops=${OPTARG};;
        i) dataInputFile=${OPTARG};;
        r) dataInputRun=${OPTARG};;
        m) maxGens=${OPTARG};;
    esac
done

# Add max gens, data input file and run to flags
flags="${flags} -M${maxGens} -I${dataInputFile} -r${dataInputRun}";

if [ "${printOutput}" = "false" ];
then
    consoleOutput="/dev/null";
fi

echo -e "[INFO] \tStarting experiment...";
echo -e "[INFO] \tFlags are ${flags}, folderBaseName: ${folderBaseName}";

### Method 1 ### -> Initialize and represent c_variables as x \in [0, 1]. Base boundaries on this representation.
# Init method 1: linear scaling of  #bins in initial population. Boundaries are based on min and max values in the data.
for ((i=4; i<= 9; i++)); do timeout 60m build/MixedIntegerGOMEA_O -P1004 -n1000 -O${folderBaseName}/${dataInputFile}/Method_1_init_1/maxDiscretizations_${i} ${flags} -X${i} >${consoleOutput}; done 
# # Init method 2: initialize all c_variables between 0 and 0.99, then normalize. (Also normalize whenever c_vars are updated)
# #           Boundaries are based on min and max values in the data.
# for ((i=2; i<= 9; i++)); do timeout 60m build/MixedIntegerGOMEA_O -P1004 -n1000 -O${folderBaseName}/${dataInputFile}/Method_1_init_2/maxDiscretizations_${i} ${flags} -X${i} -u >${consoleOutput}; done 
# Init method 3: Peter's method: force an equal amount of bins per continuous node in the population, initialize first n (for n bins) between 0 and 0.99 and then normalize.
#           The other c_vars are randomly initialized (but of course not used for boundaries). Boundaries are based on how many samples are in each bin, set halfway between last and first sample in next bin.
for ((i=4; i<= 9; i++)); do timeout 60m build/MixedIntegerGOMEA_O -P1004 -n1000 -O${folderBaseName}/${dataInputFile}/Method_1_init_3/maxDiscretizations_${i} ${flags} -X${i} -g >${consoleOutput}; done 

wait;
echo -e "[RUNNING] \tMethod 1 experiments done, starting method 3 and 3* experiments...";

### Method 3 ### -> Extra c-var per continuous node, representing the number of bins. 
# Init method 1: linear scaling of  #bins in initial population. Boundaries are based on min and max values in the data.
for ((i=4; i<= 9; i++)); do timeout 60m build/MixedIntegerGOMEA_O -P1004 -n1000 -O${folderBaseName}/${dataInputFile}/Method_3_init_1/maxDiscretizations_${i} ${flags} -X${i} -Y >${consoleOutput}; done 
# # Init method 2: initialize all c_variables between 0 and 0.99, then normalize. (Also normalize whenever c_vars are updated)
# #           Boundaries are based on min and max values in the data.
# for ((i=2; i<= 9; i++)); do timeout 60m build/MixedIntegerGOMEA_O -P1004 -n1000 -O${folderBaseName}/${dataInputFile}/Method_1_init_2/maxDiscretizations_${i} ${flags} -X${i} -u >${consoleOutput}; done 
# Init method 3: Peter's method: force an equal amount of bins per continuous node in the population, initialize first n (for n bins) between 0 and 0.99 and then normalize.
#           The other c_vars are randomly initialized (but of course not used for boundaries). Boundaries are based on how many samples are in each bin, set halfway between last and first sample in next bin.
for ((i=4; i<= 9; i++)); do timeout 60m build/MixedIntegerGOMEA_O -P1004 -n1000 -O${folderBaseName}/${dataInputFile}/Method_3_init_3/maxDiscretizations_${i} ${flags} -X${i} -g -Y >${consoleOutput}; done 

### Method 3* ### -> Extra c-var per continuous node, representing the number of bins. Added flag to force all the assigned bins to be used. 
# Init method 1: linear scaling of  #bins in initial population. Boundaries are based on min and max values in the data.
for ((i=4; i<= 9; i++)); do timeout 60m build/MixedIntegerGOMEA_O -P1004 -n1000 -O${folderBaseName}/${dataInputFile}/Method_3star_init_1/maxDiscretizations_${i} ${flags} -X${i} -Y -f >${consoleOutput}; done 
# # Init method 2: initialize all c_variables between 0 and 0.99, then normalize. (Also normalize whenever c_vars are updated)
# #           Boundaries are based on min and max values in the data.
# for ((i=2; i<= 9; i++)); do timeout 60m build/MixedIntegerGOMEA_O -P1004 -n1000 -O${folderBaseName}/${dataInputFile}/Method_1_init_2/maxDiscretizations_${i} ${flags} -X${i} -u >${consoleOutput}; done 
# Init method 3: Peter's method: force an equal amount of bins per continuous node in the population, initialize first n (for n bins) between 0 and 0.99 and then normalize.
#           The other c_vars are randomly initialized (but of course not used for boundaries). Boundaries are based on how many samples are in each bin, set halfway between last and first sample in next bin.
for ((i=4; i<= 9; i++)); do timeout 60m build/MixedIntegerGOMEA_O -P1004 -n1000 -O${folderBaseName}/${dataInputFile}/Method_3star_init_3/maxDiscretizations_${i} ${flags} -X${i} -g -Y -f >${consoleOutput}; done 
wait;
echo -e "[RUNNING] \tDone running experiments!";