#!/bin/bash

lmc=9.5
zeta=0.8
term=earth

PAR_DIR=/u/kgrunthal/HD/par/isotropic/
MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_RNCGW$lmc\_zeta$zeta\/
#MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_RNCGW$lmc\_$term\/

#MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_GWB_test/

RESULT_DIR=/u/kgrunthal/HD/out/

OUTFILE=OS_spectrum_RNCGW$lmc\_zeta$zeta
#OUTFILE=OS_spectrum_RNCGW$lmc\_$term

#NAME=RNCGW$lmc\_zeta$zeta
NAME=RNCGW$lmc\_$term

if [ ! -d "$MCMC_OUTDIR/" ]; then
  mkdir $MCMC_OUTDIR
fi

if [ ! -d "$MCMC_OUTDIR/slurm_output" ]; then
  mkdir $MCMC_OUTDIR/slurm_output
fi


#OUTFILE=testing

#sbatch -p short.q --time=04:00:00 --mem=8GB --output=$MCMC_OUTDIR/slurm_output/stdout_analysis.out --error=$MCMC_OUTDIR/slurm_output/stderr_analysis.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/EPTA_singularity/ python3 /u/kgrunthal/HD/HD_with_OS_analysis.py --par $PAR_DIR --outdir $MCMC_OUTDIR --result $RESULT_DIR/$OUTFILE --psrpickle $MCMC_OUTDIR/psrs.pkl --signal RN,CRN_fs --ptamodel TM,WN,RN,CRN_fs --lmc 9.0 --fgw 22.3 --ncgw 1 --zeta 1.0 --OStype spectrum --N 1000 --psrTerm"

sbatch -p short.q --time=04:00:00 --mem=8GB --output=$MCMC_OUTDIR/slurm_output/stdout_analysis.out --error=$MCMC_OUTDIR/slurm_output/stderr_analysis.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/EPTA_singularity/ python3 /u/kgrunthal/HD/HD_with_OS_analysis.py --par $PAR_DIR --outdir $MCMC_OUTDIR --result $RESULT_DIR/$OUTFILE --psrpickle $MCMC_OUTDIR/psrs.pkl --ptamodel TM,WN,RN,CRN_fs --OStype spectrum --N 1000 --psrTerm"

#singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/EPTA_singularity/ python3 /u/kgrunthal/HD/HD_with_OS_analysis.py --par $PAR_DIR --outdir $MCMC_OUTDIR/powerlaw/ --result $RESULT_DIR/$OUTFILE --psrpickle $MCMC_OUTDIR/psrs.pkl --signal CRN_pl --ptamodel TM,WN,CRN_pl --lmc 9.0 --fgw 22.3 --ncgw 1 --zeta 1.0 --OStype powerlaw --N 1000 --psrTerm
