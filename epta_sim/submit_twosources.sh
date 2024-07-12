#!/bin/bash

sbatch -o=./test_GW_twobins_7-10/out.txt -e=./test_GW_twobins_7-10/err.txt -q short.q --time=00:30:00 --mem=12GB --wrap="singularity exec -B /scratch/,/u/,/hercules/ /scratch/kgrunthal/EPTA_singularity/ python3 ../libstempo_multibin.py"


