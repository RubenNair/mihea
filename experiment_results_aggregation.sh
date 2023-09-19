#!/bin/bash
RED=$(tput setaf 1)
GREEN=$(tput setaf 2)
NC='\e[0m' # No Color
folder="./output/test/GAMBIT_K/popsizeMult_1.00000/2_runs"
function process () {
    count=0;
    total=0;
    totalEvals=0;
    allEvals=();
    popsize=-1;
    {
        read
        while IFS= read -r line; do
        if [[ ${line##* } -eq "1" ]]; then
            count=$((count+1));
            evals=${line%% *};
            totalEvals=$((totalEvals+evals));
            popsize=$(echo $line | awk -F" " '{print $(NF-1)}')
            allEvals+=($evals)
        fi
        total=$((total+1))
        done 
    } < "$1"
    if [[ $count -gt $((total - 2)) ]]; then
        avgEvals=$(bc <<< "scale=2; $totalEvals/$count")
        echo -e "${GREEN}$1: ${count}/${total}${NC} popsize: ${popsize} avgEvals: ${avgEvals} allEvals: ${allEvals[@]}}"
    else
        echo -e "${RED}$1: ${count}/${total}${NC} popsize: ${popsize}"
    fi
}
export -f process

while getopts f: flag
do
    case "${flag}" in
        f) folder=${OPTARG};;
    esac
done

find $folder -name \*statistics.txt -exec bash -c 'process "$@"' bash {} \;

