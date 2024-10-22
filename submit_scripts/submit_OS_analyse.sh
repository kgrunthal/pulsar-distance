#!/bin/bash

PAR_DIR=/u/kgrunthal/HD/par/isotropic/
RESULT_DIR=/u/kgrunthal/HD/out/


: <<'END_COMMENT'
for i in {1..10}; do
for lmc in 8.5 9.0 9.5 ; do
    for zeta in 0.8 0.9 1.0 ; do
        RESULT_DIR=/u/kgrunthal/HD/out/WN_CGW/run_$i\/
        MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_CGW$lmc\_zeta$zeta\_$i\/
        #MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_RNCGW$lmc\_earth/
        OUTFILE=OS_spectrum_CGW$lmc\_zeta$zeta
        NAME=$lmc\_$zeta\_spec

        if [ ! -d "$MCMC_OUTDIR/" ]; then
            mkdir $MCMC_OUTDIR
        fi

        if [ ! -d "$MCMC_OUTDIR/slurm_output" ]; then
            mkdir $MCMC_OUTDIR/slurm_output
        fi

        sbatch -p short.q --time=04:00:00 --mem=8GB --output=$MCMC_OUTDIR/slurm_output/stdout_analysis.out --error=$MCMC_OUTDIR/slurm_output/stderr_analysis.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/EPTA_singularity/ python3 /u/kgrunthal/HD/HD_with_OS_analysis.py --par $PAR_DIR --outdir $MCMC_OUTDIR --result $RESULT_DIR/$OUTFILE --psrpickle $MCMC_OUTDIR/psrs.pkl --ptamodel TM,WN,RN,CRN_fs --OStype spectrum --N 1000 --psrTerm"

    done
done
done
END_COMMENT


for i in {11..50}; do
for lmc in 8.5 9.0 9.5 ; do
    for pd in 1.0 1.5 2.0 ; do
        RESULT_DIR=/u/kgrunthal/HD/out/20PSR-galactic/
        MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_IPTA20_CGW$lmc\_pd$pd\_$i\/
        #MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_CGW$lmc\_earth/
        OUTFILE=OS_spectrum_IPTA20_CGW$lmc\_pd$pd
        NAME=$lmc\_$pd\_spec

        if [ ! -d "$MCMC_OUTDIR/" ]; then
            mkdir $MCMC_OUTDIR
        fi

        if [ ! -d "$MCMC_OUTDIR/slurm_output" ]; then
            mkdir $MCMC_OUTDIR/slurm_output
        fi

        if [ ! -d "$RESULT_DIR/" ]; then
            mkdir $RESULT_DIR
        fi

        if [ ! -d "$RESULT_DIR/run_$i/" ]; then
            mkdir $RESULT_DIR/run_$i/
        fi

        sbatch -p long.q --time=17:00:00 --mem=12GB --output=$MCMC_OUTDIR/slurm_output/analysis.out --error=$MCMC_OUTDIR/slurm_output/analysis.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/MPTA_singularity/ python3 /u/kgrunthal/HD/HD_with_OS_analysis.py --par $PAR_DIR --outdir $MCMC_OUTDIR/ --result $RESULT_DIR/run_$i\/$OUTFILE --psrpickle $MCMC_OUTDIR/psrs.pkl --ptamodel TM,WN,CRN_fs --OStype spectrum --N 1000"

    done
done
done



#sbatch -p short.q --time=04:00:00 --mem=8GB --output=$MCMC_OUTDIR/slurm_output/stdout_analysis.out --error=$MCMC_OUTDIR/slurm_output/stderr_analysis.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/EPTA_singularity/ python3 /u/kgrunthal/HD/HD_with_OS_analysis.py --par $PAR_DIR --outdir $MCMC_OUTDIR --result $RESULT_DIR/$OUTFILE --psrpickle $MCMC_OUTDIR/psrs.pkl --signal RN,CRN_fs --ptamodel TM,WN,RN,CRN_fs --lmc 9.0 --fgw 22.3 --ncgw 1 --zeta 1.0 --OStype spectrum --N 1000 --psrTerm"



#term=earth


#OUTFILE=OS_spectrum_RNCGW$lmc\_$term
#NAME=RNCGW$lmc\_$term
#MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_RNCGW$lmc\_$term\/

