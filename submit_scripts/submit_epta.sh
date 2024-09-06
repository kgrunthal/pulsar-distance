#!/bin/bash



sbatch -p short.q --time=01:00:00 --mem=12GB --output=./epta_sim/slurm_output/spectrum.out --error=./epta_sim/slurm_output/spectrum.err --job-name=EPTA --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /u/kgrunthal/EPTA_ENTERPRISE.sif python3 /u/kgrunthal/HD/epta_analysis.py"
