#!/bin/bash


for i in {1..5}; do
GWBtype=GWBsinglebin_lastbin
run=$i

PAR_DIR=/u/kgrunthal/HD/par/isotropic/
MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_$GWBtype\_$run\/
RESULT_DIR=/u/kgrunthal/HD/out/
NAME=HD_GWB


if [ ! -d "$MCMC_OUTDIR/freespectrum/" ]; then
  mkdir $MCMC_OUTDIR/freespectrum/
  mkdir $MCMC_OUTDIR/powerlaw/
fi

if [ ! -d "$MCMC_OUTDIR/slurm_output" ]; then
  mkdir $MCMC_OUTDIR/slurm_output
fi


OUTFILE=OS_spectrum_$GWBtype\_$run

sbatch -p long.q --time=16:00:00 --mem=8GB --output=$MCMC_OUTDIR/slurm_output/spectrum.out --error=$MCMC_OUTDIR/slurm_output/spectrum.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /u/kgrunthal/EPTA_ENTERPRISE.sif python3 /u/kgrunthal/HD/HD_with_OS.py --par $PAR_DIR --outdir $MCMC_OUTDIR/freespectrum/ --result $RESULT_DIR/$OUTFILE --signal GWB,CRN --lmc 9.0 --fgw 22.3 --ncgw 1 --zeta 0.8 --psrTerm --psrpickle $MCMC_OUTDIR/psrs.pkl --ptamodel TM,WN,CRN_fs --OStype spectrum --N 1000"


OUTFILE=OS_powerlaw_$GWBtype\_$run

sbatch -p long.q --time=16:00:00 --mem=8GB --output=$MCMC_OUTDIR/slurm_output/powerlaw.out --error=$MCMC_OUTDIR/slurm_output/powerlaw.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /u/kgrunthal/EPTA_ENTERPRISE.sif python3 /u/kgrunthal/HD/HD_with_OS.py --par $PAR_DIR --outdir $MCMC_OUTDIR/powerlaw --result $RESULT_DIR/$OUTFILE --signal GWB,CRN --lmc 9.0 --fgw 22.3 --ncgw 1 --zeta 0.8 --psrTerm --psrpickle $MCMC_OUTDIR/psrs.pkl --ptamodel TM,WN,CRN_pl --OStype powerlaw --N 1000"
done
