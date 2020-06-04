#!/bin/bash

UsedCores=3

export LD_LIBRARY_PATH=/opt/OpenBLAS/lib/
export OPENBLAS_NUM_THREADS=1


array=(`find . -name "*run"`);
for run in "${array[@]}"
do    
    echo "start --> $(dirname $run)";
    cd $(dirname $run) && nohup ./run > screen.out && cd .. && cd .. &

    background=( $(jobs -p) )
    if (( ${#background[@]} == UsedCores )); then
        wait -n
    fi  
done

