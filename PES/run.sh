#!/bin/bash

cores=3

export OPENBLAS_NUM_THREADS=1

array=(`find . -name "*run"`);
for run in "${array[@]}"
do    
    echo run
    cd $(dirname $run) && nohup ./run > screen.out && cd .. && cd .. &

    background=( $(jobs -p) )
    if (( ${#background[@]} == cores )); then
        wait -n
    fi  
done

