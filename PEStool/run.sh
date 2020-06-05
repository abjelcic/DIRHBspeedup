#!/bin/bash

NoSimultaneousTasks=16

export LD_LIBRARY_PATH=/opt/OpenBLAS/lib/
export OPENBLAS_NUM_THREADS=2


array=(`find . -name "*run"`);
for run in "${array[@]}"
do    
    echo "start --> $(dirname $run)";
    cd $(dirname $run) && nohup ./run > screen.out && cd .. && cd .. &

    background=( $(jobs -p) )
    if (( ${#background[@]} == NoSimultaneousTasks )); then
        wait -n
    fi  
done
