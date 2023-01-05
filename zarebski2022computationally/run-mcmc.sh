#!/usr/bin/env bash


chunk_size=10                   # the number of processes in each chunk
limit=51                        # the total number of processes (plus one)
ix=1                            # the number to start on

outer_limit="$(($limit-$chunk_size))"
while [ $ix -le $outer_limit ]
do
    pids=()
    chunk_limit="$(($ix+$chunk_size))"
    while [ $ix -lt $chunk_limit ]
    do
	      padded_num=$(printf "%03d" $ix) # because we want zero padded numbers
        ant -DbeastXML=out/s1/timtam-scenario-1-sample-$padded_num.xml mcmc & pids+=($!)
	      ((ix++))
    done
    wait "${pids[@]}"
    echo "---"
done
