#!/bin/bash

PAR_DIR=/u/kgrunthal/HD/par/isotropic/
MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_GWB_test/
RESULT_DIR=/u/kgrunthal/HD/out/
NAME=HD_GWB_new


if [ ! -d "$MCMC_OUTDIR/freespectrum/" ]; then
  mkdir $MCMC_OUTDIR/freespectrum/
  mkdir $MCMC_OUTDIR/powerlaw/
fi

if [ ! -d "$MCMC_OUTDIR/slurm_output" ]; then
  mkdir $MCMC_OUTDIR/slurm_output
fi


OUTFILE=OS_spectrum_GWB_test

sbatch -p short.q --time=04:00:00 --mem=8GB --output=$MCMC_OUTDIR/slurm_output/spectrum.out --error=$MCMC_OUTDIR/slurm_output/spectrum.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /u/kgrunthal/EPTA_ENTERPRISE.sif python3 /u/kgrunthal/HD/HD_with_OS_analysis.py --par $PAR_DIR --outdir $MCMC_OUTDIR/freespectrum/ --result $RESULT_DIR/$OUTFILE --signal GWB,CRN --lmc 9.0 --fgw 22.3 --ncgw 1 --zeta 0.8 --psrTerm --psrpickle $MCMC_OUTDIR/psrs.pkl --N 1000 --ptamodel TM,WN,CRN_fs --OStype spectrum"


OUTFILE=OS_powerlaw_GWB_test

sbatch -p short.q --time=00:30:00 --mem=8GB --output=$MCMC_OUTDIR/slurm_output/powerlaw.out --error=$MCMC_OUTDIR/slurm_output/powerlaw.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /u/kgrunthal/EPTA_ENTERPRISE.sif python3 /u/kgrunthal/HD/HD_with_OS_analysis.py --par $PAR_DIR --outdir $MCMC_OUTDIR/powerlaw --result $RESULT_DIR/$OUTFILE --signal GWB,CRN --lmc 9.0 --fgw 22.3 --ncgw 1 --zeta 0.8 --psrTerm --psrpickle $MCMC_OUTDIR/psrs.pkl --N 1000 --ptamodel TM,WN,CRN_pl --OStype powerlaw"

