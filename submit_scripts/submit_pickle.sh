#!/bin/bash

#PAR_DIR=/u/kgrunthal/HD/par/ring_20/
#PAR_DIR=/u/kgrunthal/HD/par/isotropic/
#PAR_DIR=/u/kgrunthal/HD/ska_sim/par_20/
PAR_DIR=/u/kgrunthal/HD/ska_sim/par_4/

OUTFILE=deleteme.txt
RESULT_DIR=/u/kgrunthal/HD/out/

: <<'END_COMMENT'
for i in {51..100} ; do
for lmc in 8.5 9.0 9.5 ; do
    for zeta in 0.8 0.9 1.0 ; do
        MCMC_OUTDIR=/u/kgrunthal/HD/OS_simulations/MCMCout_ring_CGW$lmc\_zeta$zeta\_$i\/
        #MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_RNCGW$lmc\_earth/
        NAME=$lmc\_$zeta

        if [ ! -d "$MCMC_OUTDIR/" ]; then
            mkdir $MCMC_OUTDIR
        fi

        if [ ! -d "$MCMC_OUTDIR/slurm_output" ]; then
            mkdir $MCMC_OUTDIR/slurm_output
        fi

        sbatch -p short.q --time=00:15:00 --mem=8GB --output=$MCMC_OUTDIR/slurm_output/pickle.out --error=$MCMC_OUTDIR/slurm_output/pickle.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/EPTA_singularity/ python3 /u/kgrunthal/HD/make_pickle.py --par $PAR_DIR --outdir $MCMC_OUTDIR --result $RESULT_DIR/$OUTFILE --fgw 22.3 --ncgw 1 --signal CGW --lmc $lmc --zeta $zeta --psrTerm"

    done
done
done
END_COMMENT

#####################

#: <<'END_COMMENT'
for i in {1..50} ; do
#for lmc in 8.5 9.0; do
#    for pd in 1.0 1.5 2.0; do
for lmc in 8.7; do    
    for pd in 1.0 4.0; do
        #MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_isotropic_earth_CGW$lmc\_pd$pd\_$i\/
        #MCMC_OUTDIR=/u/kgrunthal/HD/OS_simulations/isotropic_20/MCMCout_galactic_CGW$lmc\_pd$pd\_$i/
        #MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_galactic20_300ns_CGW$lmc\_pd$pd\_$i\/
        #MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_isotropic_scaleto8.7_CGW$lmc\_pd$pd\_$i/
        #MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_isotropic_5h_CGW$lmc\_pd$pd\_$i/
        MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_galactic4_scaleto9.0_CGW$lmc\_pd$pd\_$i/
        NAME=pickle_$lmc\_$pd

        if [ ! -d "$MCMC_OUTDIR/" ]; then
            mkdir $MCMC_OUTDIR
        fi

        if [ ! -d "$MCMC_OUTDIR/slurm_output" ]; then
            mkdir $MCMC_OUTDIR/slurm_output
        fi

        # normal
        #sbatch -p short.q --time=00:15:00 --mem=8GB --output=$MCMC_OUTDIR/slurm_output/pickle.out --error=$MCMC_OUTDIR/slurm_output/pickle.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/EPTA_singularity/ python3 /u/kgrunthal/HD/make_pickle.py --par $PAR_DIR --outdir $MCMC_OUTDIR --result $RESULT_DIR/$OUTFILE --fgw 22.3 --ncgw 1 --signal CGW --lmc $lmc --zeta -1 --pdistance $pd --pdist_fix --psrTerm"

        # scale to 8.7
        sbatch -p short.q --time=00:15:00 --mem=8GB --output=$MCMC_OUTDIR/slurm_output/pickle.out --error=$MCMC_OUTDIR/slurm_output/pickle.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/EPTA_singularity/ python3 /u/kgrunthal/HD/make_pickle.py --par $PAR_DIR --outdir $MCMC_OUTDIR --result $RESULT_DIR/$OUTFILE --fgw 22.3 --ncgw 1 --signal CGW --lmc $lmc --zeta -1 --pdistance $pd --pdist_fix --psrTerm --scale_to 9.0"

        # 5*h of 8.7 at 15 kpc
        #sbatch -p short.q --time=00:15:00 --mem=8GB --output=$MCMC_OUTDIR/slurm_output/pickle.out --error=$MCMC_OUTDIR/slurm_output/pickle.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/EPTA_singularity/ python3 /u/kgrunthal/HD/make_pickle.py --par $PAR_DIR --outdir $MCMC_OUTDIR --result $RESULT_DIR/$OUTFILE --fgw 22.3 --ncgw 1 --signal CGW --lmc $lmc --zeta -1 --pdistance $pd --pdist_fix --psrTerm --scale_to 9.119"
    done
done
done
#END_COMMENT

######################
PAR_DIR=/u/kgrunthal/HD/par/isotropic/
: <<'END_COMMENT'
for i in {1..5} ;do
MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_GWBsinglebin_lastbin_$i\/
NAME=pickle

if [ ! -d "$MCMC_OUTDIR/" ]; then
    mkdir $MCMC_OUTDIR
fi

if [ ! -d "$MCMC_OUTDIR/slurm_output" ]; then
    mkdir $MCMC_OUTDIR/slurm_output
fi


sbatch -p short.q --time=00:15:00 --mem=8GB --output=$MCMC_OUTDIR/slurm_output/pickle.out --error=$MCMC_OUTDIR/slurm_output/pickle.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/EPTA_singularity/ python3 /u/kgrunthal/HD/make_pickle.py --par $PAR_DIR --outdir $MCMC_OUTDIR --result $RESULT_DIR/$OUTFILE --signal GWBsinglebin"
done
END_COMMENT

