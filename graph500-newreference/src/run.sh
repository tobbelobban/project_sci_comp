#!/bin/bash
SCALE=8
EFACTOR=10
export TMPFILE=${SCALE}_${EFACTOR}.bin
export REUSEFILE=1
mpirun ./graph500_reference_bfs ${SCALE} ${EFACTOR}
