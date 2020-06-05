#!/bin/bash

NoSimultaneousTasks=16
export OPENBLAS_NUM_THREADS=2


export LD_LIBRARY_PATH=/opt/OpenBLAS/lib/

array=(`find . -name "*.pes"`);
for run in "${array[@]}"
do    
    echo "start --> $(dirname $run)";
    cd $(dirname $run) && nohup "./$(basename $run)" > screen.out && cd .. && cd .. &

    background=( $(jobs -p) )
    if (( ${#background[@]} == NoSimultaneousTasks )); then
        wait -n
    fi  
done

