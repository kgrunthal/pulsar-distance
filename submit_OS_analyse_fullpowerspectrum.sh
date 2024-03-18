#!/bin/bash

PAR_DIR=/u/kgrunthal/HD/par/isotropic/
MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_GWBbroken/
RESULT_DIR=/u/kgrunthal/HD/out/
OUTFILE=OS_spectrum_GWBbroken_fullspectrum.txt
NAME=GWB


#sbatch -p short.q --time=04:00:00 --mem=8GB --output=$MCMC_OUTDIR/slurm_output/stdout_analysis.out --error=$MCMC_OUTDIR/slurm_output/stderr_analysis.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/EPTA_singularity/ python3 /u/kgrunthal/HD/HD_with_OS_fullpowerspectrum.py --par $PAR_DIR --outdir $MCMC_OUTDIR --result $RESULT_DIR/$OUTFILE --signal GWB --ptamodel TM,WN,CRN --lmc 9.0 --fgw 22.3 --ncgw 1 --zeta 0.8 --noisemarginalised --N 1000 --psrTerm --psrpickle $MCMC_OUTDIR/psrs.pkl"


sbatch -p short.q --time=00:10:00 --mem=8GB --output=$MCMC_OUTDIR/slurm_output/fullspectrum.out --error=$MCMC_OUTDIR/slurm_output/fullspectrum.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/EPTA_singularity/ python3 /u/kgrunthal/HD/HD_with_OS_fullpowerspectrum.py --par $PAR_DIR --outdir $MCMC_OUTDIR/freespectrum/ --result $RESULT_DIR/$OUTFILE --signal GWB --ptamodel TM,WN,CRN --lmc 9.0 --fgw 22.3 --ncgw 1 --zeta 0.8 --N 1000 --psrTerm --psrpickle $MCMC_OUTDIR/psrs.pkl"
