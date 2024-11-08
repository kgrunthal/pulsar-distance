#!/bin/bash

#:<< END_COMMENT
for i in {1..100}; do
for lmc in 9.5; do
#for lmc in 8.7; do
for pd in 1.0 2.0; do
#for pd in 1.0 4.0; do
    #MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_isotropic_newparams_CGW$lmc\_pd$pd\_$i\/
    #MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_isotropic_fixpd_CGW$lmc\_pd$pd\_$i\/
    MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_galactic20_CGW$lmc\_pd$pd\_$i\/

    NAME=CGW_$lmc\_$pd

    RESULT_DIR_noPT=$MCMC_OUTDIR/CGWsearch_noPT/
    RESULT_DIR_noPT_fixparam=$MCMC_OUTDIR/CGWsearch_noPT_fixPhi0Psi/


    if [ ! -d "$RESULT_DIR_noPT/" ]; then
         mkdir $RESULT_DIR_noPT
    fi

    #if [ ! -d "$RESULT_DIR_noPT/" ]; then
    #     mkdir $RESULT_DIR_noPT_fixPhi0Psi/
    #fi


    ### all parameters ###
    # sampling
    sbatch -p short.q --time=04:00:00 --mem=12GB --output=$MCMC_OUTDIR/slurm_output/CGW_noPT.out --error=$MCMC_OUTDIR/slurm_output/CGW_noPT.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/MPTA_singularity/ python3 CGW_search.py --basedir $MCMC_OUTDIR --outdir $RESULT_DIR_noPT --ptamodel TM,WN,CGW --analysis --use_distance --run_sampler"

    # analysis only 
#    sbatch -p short.q --time=00:10:00 --mem=12GB --output=$MCMC_OUTDIR/slurm_output/chain_analysis.out --error=$MCMC_OUTDIR/slurm_output/chain_analysis.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/MPTA_singularity/ python3 CGW_search.py --basedir $MCMC_OUTDIR --outdir $RESULT_DIR_noPT --ptamodel TM,WN,CGW --analysis --use_distance"


    ### fix Phi0 and Psi ### 
    # sampling
#    sbatch -p short.q --time=04:00:00 --mem=12GB --output=$MCMC_OUTDIR/slurm_output/CGW_noPT_fixPhi0Psi.out --error=$MCMC_OUTDIR/slurm_output/CGW_noPT_fixPhi0Psi.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/MPTA_singularity/ python3 CGW_search.py --basedir $MCMC_OUTDIR --outdir $RESULT_DIR_noPT --ptamodel TM,WN,CGW --analysis --use_distance --run_sampler"

done
done
done
#END_COMMENT



:<< END_COMMENT
# single launches
lmc=9.5
pd=1.5
i=38

MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_CGW$lmc\_pd$pd\_$i\/
NAME=CGW_$lmc\_$pd

RESULT_DIR_noPT=$MCMC_OUTDIR/CGWsearch_noPT/
RESULT_DIR_PT=$MCMC_OUTDIR/CGWsearch_PT/

if [ ! -d "$RESULT_DIR_noPT/" ]; then
    mkdir $RESULT_DIR_noPT
fi

if [ ! -d "$RESULT_DIR_PT/" ]; then
    mkdir $RESULT_DIR_PT
fi

#sbatch -p short.q --time=02:00:00 --mem=12GB --output=$MCMC_OUTDIR/slurm_output/CGW_noPT.out --error=$MCMC_OUTDIR/slurm_output/CGW_noPT.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /u/kgrunthal/EPTA_ENTERPRISE.sif python3 CGW_search.py --basedir $MCMC_OUTDIR --outdir $RESULT_DIR_noPT --ptamodel TM,WN,CGW --analysis"

sbatch -p short.q --time=02:00:00 --mem=12GB --output=$MCMC_OUTDIR/slurm_output/CGW_PT.out --error=$MCMC_OUTDIR/slurm_output/CGW_PT.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /u/kgrunthal/EPTA_ENTERPRISE.sif python3 CGW_search.py --basedir $MCMC_OUTDIR --outdir $RESULT_DIR_PT --ptamodel TM,WN,CGW --psrTerm --analysis"
END_COMMENT

