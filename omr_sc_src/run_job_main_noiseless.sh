#!/bin/bash

# if you have access to a comptuer cluster, please modify the script and use "qsub" to submit the jobs to the cluster so that they can be ran in parallel
mkdir -p ./results_noiseless

for l in 10 ; do
    for p in 100 ; do
        for r in 3; do
            for m in $(seq 1 10); do
                ./run_main_noiseless.sh ${l} ${p} ${r} ${m} 
            done
        done
    done
done
