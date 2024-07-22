#!/bin/bash

FOLDER=test_GW7_RN12
sbatch -o=./$FOLDER\/out.txt -e=./$FOLDER\/err.txt -p short.q --time=00:30:00 --mem=12GB --wrap="singularity exec -B /scratch/,/u/,/hercules/ /scratch/kgrunthal/EPTA_singularity/ python3 ../libstempo_multibin.py"


