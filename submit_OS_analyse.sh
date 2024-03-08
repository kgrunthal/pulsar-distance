#!/bin/bash

PAR_DIR=/u/kgrunthal/HD/par/isotropic/
MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_GWB/
RESULT_DIR=/u/kgrunthal/HD/out/noisemarginalised/
OUTFILE=OS_spectrum_WNGWB_NM.txt
NAME=GWB

mkdir $MCMC_OUTDIR

sbatch -p short.q --time=04:00:00 --mem=8GB --output=$MCMC_OUTDIR/slurm_output/stdout_analysis.out --error=$MCMC_OUTDIR/slurm_output/stderr_analysis.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/EPTA_singularity/ python3 /u/kgrunthal/HD/HD_with_OS_analysis.py --par $PAR_DIR --outdir $MCMC_OUTDIR --result $RESULT_DIR/$OUTFILE --signal GWB --ptamodel TM,WN,CRN --lmc 9.0 --fgw 22.3 --ncgw 1 --zeta 0.8 --noisemarginalised --N 10000 --psrTerm"
