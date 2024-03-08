#!/bin/bash

PAR_DIR=/u/kgrunthal/HD/par/isotropic/
MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_RNCGW9.0_zeta1.0/
RESULT_DIR=/u/kgrunthal/HD/out/
OUTFILE=OS_spectrum_GWB_CRN.txt
NAME=HD_RNCGW


mkdir $MCMC_OUTDIR
mkdir $MCMC_OUTDIR/slurm_output



sbatch -p short.q --time=00:15:00 --mem=8GB --output=$MCMC_OUTDIR/slurm_output/pickle.out --error=$MCMC_OUTDIR/slurm_output/pickle.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /u/kgrunthal/EPTA_ENTERPRISE.sif python3 /u/kgrunthal/HD/make_pickle.py --par $PAR_DIR --outdir $MCMC_OUTDIR --result $RESULT_DIR/$OUTFILE --signal RN,CGW --lmc 9.0 --fgw 22.3 --ncgw 1 --zeta 1.0 --psrTerm"
