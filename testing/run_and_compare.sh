#!/bin/bash

mkdir -p simulation/output

error1=$(python3 simulation/main.py -s 888 -n 20 -t 3 -o testing/output/test1/ 2>&1 >testing/output/test1/stdout)
echo "Test 1 diff:"
diff testing/output/test1/final_state.csv testing/correct_outputs/test1_final_state.csv

cp testing/inputs/test1_snapshots.hdf5 testing/inputs/loadsnap_input.hdf5
error_loadsnap=$(python3 simulation/main.py -s 888 -n 20 -T 2.0 -t 3 -o testing/output/loadsnap/ -i testing/inputs 2>&1 >testing/output/loadsnap/stdout)
echo "Snapshot load test diff:"
diff testing/output/loadsnap/final_state.csv testing/correct_outputs/test1_final_state.csv

echo "Errors from test 1:"
echo $error1

echo "Errors from snapshot load test:"
echo $error_loadsnap
