#!/bin/bash

mkdir -p simulation/output

error1=$(python3 simulation/main.py -s 888 -n 20 -t 3 2>&1 >simulation/output/stdout)
echo "Test 1 diff:"
diff simulation/output/final_state.csv testing/correct_outputs/test1_final_state.csv

error_loadsnap=$(python3 simulation/main.py -s 888 -n 20 -T 2 -t 3 -i testing/inputs/test1_snapshots.hdf5 2>&1 >simulation/output/stdout)
echo "Snapshot load test diff:"
diff simulation/output/final_state.csv testing/correct_outputs/test1_final_state.csv

echo "Errors from test 1:"
echo $error1

echo "Errors from snapshot load test:"
echo $error_loadsnap
