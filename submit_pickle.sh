#!/bin/bash

lmc=9.0
zeta=1.0


PAR_DIR=/u/kgrunthal/HD/par/isotropic/
#MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_RNCGW$lmc\_zeta$zeta\/
MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_RNCGW$lmc\_earth/
#MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_GWBbroken_2_new/
RESULT_DIR=/u/kgrunthal/HD/out/

OUTFILE=OS_spectrum_GWB_CRN.txt
NAME=pickle



if [ ! -d "$MCMC_OUTDIR/" ]; then
    mkdir $MCMC_OUTDIR
fi


if [ ! -d "$MCMC_OUTDIR/slurm_output" ]; then
    mkdir $MCMC_OUTDIR/slurm_output
fi


#sbatch -p short.q --time=00:15:00 --mem=8GB --output=$MCMC_OUTDIR/slurm_output/pickle.out --error=$MCMC_OUTDIR/slurm_output/pickle.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/EPTA_singularity/ python3 /u/kgrunthal/HD/make_pickle.py --par $PAR_DIR --outdir $MCMC_OUTDIR --result $RESULT_DIR/$OUTFILE --signal GWBbroken"

sbatch -p short.q --time=00:15:00 --mem=8GB --output=$MCMC_OUTDIR/slurm_output/pickle.out --error=$MCMC_OUTDIR/slurm_output/pickle.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/EPTA_singularity/ python3 /u/kgrunthal/HD/make_pickle.py --par $PAR_DIR --outdir $MCMC_OUTDIR --result $RESULT_DIR/$OUTFILE --fgw 22.3 --ncgw 1 --signal RN,CGW --lmc $lmc --zeta $zeta"
