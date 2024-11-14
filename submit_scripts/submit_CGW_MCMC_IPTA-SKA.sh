#!/bin/bash

#:<< END_COMMENT
for i in {1..50}; do

#for lmc in 8.5 9.0 9.5; do
for lmc in 9.0; do

#for pd in low; do
#for pd in low mid high; do
for pd in full over2.0; do
#for pd in full far ipta; do
    
    #MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_IPTA_WN_CGW_RA12hDEC0deg_$lmc\_$pd\_$i\/
    MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_SKA_WN_CGW_RA12hDEC0deg_$lmc\_$pd\_$i\/
    
    NAME=$lmc\_$pd
    
    RESULT_DIR_noPT=$MCMC_OUTDIR/CGWsearch_noPT/
    RESULT_DIR_PT=$MCMC_OUTDIR/CGWsearch_PT/
    RESULT_DIR_PT_pphase=$MCMC_OUTDIR/CGWsearch_PT_pphase/

    if [ ! -d "$RESULT_DIR_noPT/" ]; then
         mkdir $RESULT_DIR_noPT
    fi

    #if [ ! -d "$RESULT_DIR_PT/" ]; then
    #     mkdir $RESULT_DIR_PT
    #fi

    #if [ ! -d "$RESULT_DIR_PT_pphase/" ]; then
    #     mkdir $RESULT_DIR_PT_pphase
    #fi

    if [ ! -f "$RESULT_DIR_noPT/maxlike.json" ]; then
        # IPTA only
        #sbatch -p long.q --time=06:00:00 --mem=20GB --output=$MCMC_OUTDIR/slurm_output/CGW_noPT.out --error=$MCMC_OUTDIR/slurm_output/CGW_noPT.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/MPTA_singularity/ python3 CGW_search_IPTA-SKA.py --basedir $MCMC_OUTDIR --outdir $RESULT_DIR_noPT --ptamodel TM,WN,CGW --noisefile /u/kgrunthal/HD/ipta_sim/WN_dictionary.json --analysis --use_distance"

        # SKA
        sbatch -p long.q --time=10:00:00 --mem=20GB --output=$MCMC_OUTDIR/slurm_output/CGW_noPT.out --error=$MCMC_OUTDIR/slurm_output/CGW_noPT.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/MPTA_singularity/ python3 CGW_search_IPTA-SKA.py --basedir $MCMC_OUTDIR --outdir $RESULT_DIR_noPT --ptamodel TM,WN,CGW --noisefile /u/kgrunthal/HD/ska_sim/FULL_WN_dictionary.json --analysis --use_distance"
    fi
 

    
#    sbatch -p short.q --time=03:00:00 --mem=12GB --output=$MCMC_OUTDIR/slurm_output/CGW_PT_pphase.out --error=$MCMC_OUTDIR/slurm_output/CGW_PT_pphase.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/MPTA_singularity/ python3 CGW_search.py --basedir $MCMC_OUTDIR --outdir $RESULT_DIR_PT_pphase --ptamodel TM,WN,CGW --pd $pd --psrTerm --analysis --use_distance" 

#    sbatch -p short.q --time=00:30:00 --mem=12GB --output=$MCMC_OUTDIR/slurm_output/CGW_PT_analysis.out --error=$MCMC_OUTDIR/slurm_output/CGW_PT_analysis.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/MPTA_singularity/ python3 CGW_search.py --basedir $MCMC_OUTDIR --outdir $RESULT_DIR_PT --ptamodel TM,WN,CGW --pd $pd --psrTerm --analysis --use_distance --sample_pdist"

#    singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/MPTA_singularity/ python3 CGW_search.py --basedir $MCMC_OUTDIR --outdir $RESULT_DIR_PT --ptamodel TM,WN,CGW --pd $pd --psrTerm --analysis --use_distance --Nsample 50000 --sample_pdist

done
done
done
#END_COMMENT



