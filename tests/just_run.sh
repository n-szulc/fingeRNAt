#!/bin/bash

mkdir -p outputs/
date > outputs/data.txt
ls -la
echo "********"
# python ../code/fingeRNAt.py -r test_inputs/3d2v.pdb -l test_inputs/redocked.sdf -f FULL -h2o -o outputs -custom test_inputs/custom-interactions.yaml -detail -debug > outputs/debug.txt
echo "********"
ls -la