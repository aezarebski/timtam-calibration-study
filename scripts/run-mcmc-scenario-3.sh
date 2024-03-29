#!/usr/bin/env bash
#
# Usage
# =====
#
# If there are 50 replicates to run for the second analysis, execute
# this script with the following command:
#
# $ bash run-mcmc-scenario-3.sh 50 2
#

# The `CHUNK_SIZE' is the variable you should edit based on your
# machine.
CHUNK_SIZE=10                   # the number of processes per chunk
NUM_REPLICATES=${1:-50}         # the total number of processes
IX=1                            # the number to start on
LIMIT=$((NUM_REPLICATES+1))

OUTER_LIMIT="$((LIMIT-CHUNK_SIZE))"

ANALYSIS_NUM=$2

CHAIN_LEN=$3

echo "Running $NUM_REPLICATES each of length $CHAIN_LEN"

while [ $IX -le $OUTER_LIMIT ]
do
    PIDS=()
    CHUNK_LIMIT=$((IX+CHUNK_SIZE))
    while [ $IX -lt $CHUNK_LIMIT ]
    do
	PADDED_NUM=$(printf "%03d" $IX) # because we want zero padded numbers
        ant -DstateFile=out/tmp/timtam-scenario-3-"$ANALYSIS_NUM"-sample-"$PADDED_NUM".xml.state -DchainLength=$CHAIN_LEN -DbeastXML=out/s3/timtam-scenario-3-"$ANALYSIS_NUM"-sample-"$PADDED_NUM".xml mcmc & PIDS+=($!)
	((IX++))
    done
    wait "${PIDS[@]}"
    echo "---"
done
