#!/bin/bash
SCALE=14
EFACTOR=20
export TMPFILE=${SCALE}_${EFACTOR}.bin
export REUSEFILE=1
mpirun ./graph500_reference_bfs ${SCALE} ${EFACTOR}
