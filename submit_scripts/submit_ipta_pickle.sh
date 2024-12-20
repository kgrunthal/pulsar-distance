#!/bin/bash

HEAD_DIR=/u/kgrunthal/HD/ipta_sim/
RESULT_DIR=/u/kgrunthal/HD/out/


for i in {1..50} ; do
for lmc in 8.7 9.0; do
#for lmc in 8.5; do    
    for pd in full over2.0; do
        MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_IPTA_WN_CGW_RA12hDEC0deg_$lmc\_$pd\_$i\/
        
        NAME=$lmc\_$pd

        if [ ! -d "$MCMC_OUTDIR/" ]; then
            mkdir $MCMC_OUTDIR
        fi

        if [ ! -d "$MCMC_OUTDIR/slurm_output" ]; then
            mkdir $MCMC_OUTDIR/slurm_output
        fi

#        sbatch -p short.q --time=00:15:00 --mem=8GB --output=$MCMC_OUTDIR/slurm_output/pickle.out --error=$MCMC_OUTDIR/slurm_output/pickle.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/EPTA_singularity/ python3 /u/kgrunthal/HD/make_pickle.py --datadir $PAR_DIR --outdir $MCMC_OUTDIR --fgw 22.3 --ncgw 1 --ptamodel CGW --lmc $lmc --zeta -1 --pdistance $pd --psrTerm"

        sbatch -p short.q --time=00:15:00 --mem=8GB --output=$MCMC_OUTDIR/slurm_output/pickle_$pd\_$i\.out --error=$MCMC_OUTDIR/slurm_output/pickle_$pd\_$i\.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/EPTA_singularity/ python3 /u/kgrunthal/HD/ipta_pickle.py --datadir $HEAD_DIR/par/all/ --outdir $MCMC_OUTDIR --ncgw 1 --zeta -1 --psr_list $HEAD_DIR/ipta_$pd\.txt --fgw 22.3 --ptamodel CGW --lmc $lmc --psrTerm"
    done
done
done
