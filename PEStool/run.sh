#!/bin/bash

NoSimultaneousTasks=4
export OPENBLAS_NUM_THREADS=1


export LD_LIBRARY_PATH=/opt/OpenBLAS/lib/




function runjob(){
    echo "start --> $(dirname $run)";
    cd $(dirname $run)
    nohup "./$(basename $run)" > screen.out
    cd ..
    cd ..
}

array=(`find . -wholename "./output/*.pes"`);
for run in "${array[@]}"
do    
    runjob &
    background=( $(jobs -p) )
    if (( ${#background[@]} == NoSimultaneousTasks )); then
        wait -n
    fi  
done

