#!/bin/bash

mkdir -p simulation/output

error1=$(python3 simulation/main.py -s 888 -n 20 -t 3 -b 0.1 2>&1 >simulation/output/stdout)
echo "Test 1 diff:"
diff simulation/output/final_state.csv testing/correct_outputs/test1_final_state.csv

error2=$(python3 simulation/main.py -s 888 -n 20 -t 3 -b 0.1 -e -1 2>&1 >simulation/output/stdout)
echo "Test 2 diff:"
diff simulation/output/final_state.csv testing/correct_outputs/test2_final_state.csv

echo "Errors (test 1 only):"
echo $error1
