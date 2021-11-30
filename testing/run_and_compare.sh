#!/bin/bash

mkdir -p simulation/output

python3 simulation/main.py -s 888 -n 20 -t 3 > simulation/output/stdout
echo "Test 1 diff:"
diff simulation/output/final_state.png testing/correct_outputs/test1_final_state.png

python3 simulation/main.py -s 888 -n 20 -t 3 -e -1 > simulation/output/stdout
echo "Test 2 diff:"
diff simulation/output/final_state.png testing/correct_outputs/test2_final_state.png
