#!/bin/bash

#PAR_DIR=/u/kgrunthal/HD/par/isotropic/
#PAR_DIR=/u/kgrunthal/HD/ipta_sim/par/all/
#PAR_DIR=/u/kgrunthal/HD/par/ring/
PAR_DIR=/u/kgrunthal/HD/ska_sim/par_20/

: <<'END_COMMENT'
# different zeta
for i in {51..100}; do
for lmc in 8.5 9.0 9.5 ; do
    for zeta in 0.8 0.9 1.0 ; do
        RESULT_DIR=/u/kgrunthal/HD/out/WN_CGW_zeta/
        MCMC_OUTDIR=/u/kgrunthal/HD/OS_simulations/ring_20/MCMCout_ring_CGW$lmc\_zeta$zeta\_$i\/
        #MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_RNCGW$lmc\_earth/
        OUTFILE=OS_spectrum_CGW$lmc\_zeta$zeta
        NAME=OS_$lmc\_$zeta

        if [ ! -d "$MCMC_OUTDIR/" ]; then
            mkdir $MCMC_OUTDIR
        fi

        if [ ! -d "$MCMC_OUTDIR/slurm_output" ]; then
            mkdir $MCMC_OUTDIR/slurm_output
        fi

        if [ ! -d "$RESULT_DIR/" ]; then
            mkdir $RESULT_DIR
        fi

        if [ ! -d "$RESULT_DIR/run_$i" ]; then
            mkdir $RESULT_DIR/run_$i
        fi
 
        sbatch -p short.q --time=02:00:00 --mem=12GB --output=$MCMC_OUTDIR/slurm_output/spectrum.out --error=$MCMC_OUTDIR/slurm_output/spectrum.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/MPTA_singularity/ python3 /u/kgrunthal/HD/HD_with_OS.py --par $PAR_DIR --outdir $MCMC_OUTDIR/ --result $RESULT_DIR/run_$i\/$OUTFILE --psrpickle $MCMC_OUTDIR/psrs.pkl --ptamodel TM,WN,CRN_fs --OStype spectrum --N 1000"

    done
done
done
END_COMMENT



################################

#: <<'END_COMMENT'
for i in {1..100}; do
#for lmc in 8.5 9.0 9.5; do
for lmc in 8.7; do
    #for pd in 1.0 1.5 2.0; do
    for pd in 1.0 4.0; do
    #for pd in full over2.0; do
        #RESULT_DIR=/u/kgrunthal/HD/out/IPTA_OS/
        #RESULT_DIR=/u/kgrunthal/HD/out/WN_CGW_pd/
        #RESULT_DIR=/u/kgrunthal/HD/out/galactic_scaleto8.7/
        RESULT_DIR=/u/kgrunthal/HD/out/galactic_h5/

        #MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_galactic20_scaleto8.7_CGW$lmc\_pd$pd\_$i\/
        MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_galactic20_5h_CGW$lmc\_pd$pd\_$i/
        #MCMC_OUTDIR=/u/kgrunthal/HD/OS_simulations/isotropic_20/MCMCout_isotropic_CGW$lmc\_pd$pd\_$i/
        #MCMC_OUTDIR=/u/kgrunthal/HD/OS_simulations/galactic_20/MCMCout_galactic_CGW$lmc\_pd$pd\_$i/
        OUTFILE=OS_spectrum_CGW$lmc\_pd$pd
        NAME=$lmc\_$pd\_spec

        if [ ! -d "$MCMC_OUTDIR/" ]; then
            mkdir $MCMC_OUTDIR
        fi

        if [ ! -d "$MCMC_OUTDIR/slurm_output" ]; then
            mkdir $MCMC_OUTDIR/slurm_output
        fi

        if [ ! -d "$MCMC_OUTDIR/PFOS" ]; then
            mkdir $MCMC_OUTDIR/PFOS
        fi

        if [ ! -d "$RESULT_DIR/" ]; then
            mkdir $RESULT_DIR
        fi

        if [ ! -d "$RESULT_DIR/run_$i" ]; then
            mkdir $RESULT_DIR/run_$i
        fi


        sbatch -p short.q --time=02:00:00 --mem=12GB --output=$MCMC_OUTDIR/slurm_output/spectrum.out --error=$MCMC_OUTDIR/slurm_output/spectrum.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /scratch/kgrunthal/MPTA_singularity/ python3 /u/kgrunthal/HD/HD_with_OS.py --par $PAR_DIR --outdir $MCMC_OUTDIR/PFOS/ --result $RESULT_DIR/run_$i\/$OUTFILE --psrpickle $MCMC_OUTDIR/psrs.pkl --ptamodel TM,WN,CRN_fs --OStype spectrum --N 1000"

    done
done
done

#END_COMMENT

#############################

: <<'END_COMMENT'

term=earth
MCMC_OUTDIR=/u/kgrunthal/HD/MCMCout_RNCGW$lmc\_$term\/
OUTFILE=OS_spectrum_RNCGW$lmc\_$term
NAME=RNCGW$lmc\_$term\/

if [ ! -d "$MCMC_OUTDIR/" ]; then
  mkdir $MCMC_OUTDIR
fi

if [ ! -d "$MCMC_OUTDIR/slurm_output" ]; then
  mkdir $MCMC_OUTDIR/slurm_output
fi



sbatch -p long.q --time=16:00:00 --mem=10GB --output=$MCMC_OUTDIR/slurm_output/spectrum.out --error=$MCMC_OUTDIR/slurm_output/spectrum.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /u/kgrunthal/EPTA_ENTERPRISE.sif python3 /u/kgrunthal/HD/HD_with_OS.py --par $PAR_DIR --outdir $MCMC_OUTDIR/ --result $RESULT_DIR/$OUTFILE --psrpickle $MCMC_OUTDIR/psrs.pkl --ptamodel TM,WN,RN,CRN_fs --OStype spectrum --N 1000"


#sbatch -p long.q --time=16:00:00 --mem=10GB --output=$MCMC_OUTDIR/slurm_output/spectrum.out --error=$MCMC_OUTDIR/slurm_output/spectrum.err --job-name=$NAME --wrap="singularity exec -B /scratch/kgrunthal/,/hercules/results/kgrunthal/,/u/kgrunthal/ /u/kgrunthal/EPTA_ENTERPRISE.sif python3 /u/kgrunthal/HD/HD_with_OS.py --par $PAR_DIR --outdir $MCMC_OUTDIR/ --result $RESULT_DIR/$OUTFILE --signal RN,CGW --lmc 9.0 --fgw 22.3 --ncgw 1 --psrTerm --psrpickle $MCMC_OUTDIR/psrs.pkl --ptamodel TM,WN,RN,CRN_fs --OStype spectrum --N 1000 --zeta 1.0 --resume"

END_COMMENT


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
