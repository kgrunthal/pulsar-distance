#!/bin/bash

FOLDER=test_bin7_A5e-12_2
sbatch -o ./$FOLDER/out.txt -e ./$FOLDER/err.txt -p short.q --time=00:30:00 --job-name=singlebin --mem=12GB --wrap="singularity exec -B /scratch/,/u/,/hercules/ /scratch/kgrunthal/EPTA_singularity/ python3 ../libstempo_singlebin.py"


