#!/bin/bash

PAR_DIR=/u/kgrunthal/HD/par/isotropic/
MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_RNCGW9.0_zeta1.0/
RESULT_DIR=/u/kgrunthal/HD/out/
NAME=RNCGW9.0_zeta1.0


if [ ! -d "$MCMC_OUTDIR/" ]; then
  mkdir $MCMC_OUTDIR
fi

if [ ! -d "$MCMC_OUTDIR/slurm_output" ]; then
  mkdir $MCMC_OUTDIR/slurm_output
fi


OUTFILE=OS_spectrum_RNCGW9.0_zeta1.0

sbatch -p long.q --time=16:00:00 --mem=8GB --output=$MCMC_OUTDIR/slurm_output/spectrum.out --error=$MCMC_OUTDIR/slurm_output/spectrum.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /u/kgrunthal/EPTA_ENTERPRISE.sif python3 /u/kgrunthal/HD/HD_with_OS.py --par $PAR_DIR --outdir $MCMC_OUTDIR/ --result $RESULT_DIR/$OUTFILE --signal RN,CGW --lmc 9.0 --fgw 22.3 --ncgw 1 --psrTerm --psrpickle $MCMC_OUTDIR/psrs.pkl --ptamodel TM,WN,RN,CRN_fs --OStype spectrum --N 1000 --zeta 1.0"




: <<'END_COMMENT'
if [ ! -d "$MCMC_OUTDIR/freespectrum/" ]; then
  mkdir $MCMC_OUTDIR/freespectrum/
  mkdir $MCMC_OUTDIR/powerlaw/
fi

if [ ! -d "$MCMC_OUTDIR/slurm_output" ]; then
  mkdir $MCMC_OUTDIR/slurm_output
fi


OUTFILE=OS_spectrum_GWBbroken

sbatch -p long.q --time=16:00:00 --mem=8GB --output=$MCMC_OUTDIR/slurm_output/spectrum.out --error=$MCMC_OUTDIR/slurm_output/spectrum.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /u/kgrunthal/EPTA_ENTERPRISE.sif python3 /u/kgrunthal/HD/HD_with_OS.py --par $PAR_DIR --outdir $MCMC_OUTDIR/freespectrum/ --result $RESULT_DIR/$OUTFILE --signal RN,CGW --lmc 9.0 --fgw 22.3 --ncgw 1 --zeta 0.8 --psrTerm --psrpickle $MCMC_OUTDIR/psrs.pkl --ptamodel TM,WN,RN,CRN_fs --OStype spectrum"


OUTFILE=OS_powerlaw_GWBbroken

sbatch -p long.q --time=16:00:00 --mem=8GB --output=$MCMC_OUTDIR/slurm_output/powerlaw.out --error=$MCMC_OUTDIR/slurm_output/powerlaw.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /u/kgrunthal/EPTA_ENTERPRISE.sif python3 /u/kgrunthal/HD/HD_with_OS.py --par $PAR_DIR --outdir $MCMC_OUTDIR/powerlaw --result $RESULT_DIR/$OUTFILE --signal GWB,CRN --lmc 9.0 --fgw 22.3 --ncgw 1 --zeta 0.8 --psrTerm --psrpickle $MCMC_OUTDIR/psrs.pkl --ptamodel TM,WN,CRN_pl --OStype powerlaw"
END_COMMENT
