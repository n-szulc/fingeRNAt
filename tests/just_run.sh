#!/bin/bash

python ../code/fingeRNAt.py -r test_inputs/3d2v.pdb -l test_inputs/redocked.sdf -f FULL -h2o -o outputs -custom test_inputs/custom-interactions.yaml -detail -debug > outputs/debug.txt