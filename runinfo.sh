#!/bin/bash

###############################################################
## Options for running PHQMD code, stabilisation & conversion #

## Steps:
export RunPhqmdCode=0           ## = 1: run PHQMD code - input: inputPHSD, output: phsd.dat (hadrons), fort.791 (cluster-baryons), fort.891 (anti-cluster-baryons)
export RunStabilisation=0       ## = 1: run stabilisation routine - input: fort.791 & fort.781, output: fort.891 & fort.881
export MakeDetectorInput=1      ## = 1: run conversion of PHQMD output into detector input - input: fort.891 & fort.881, output: unigen format

## Options for steps:
export ConvertAntiClusters=1    ## = 0: convert/stablise only clusters, = 1: convert clusters and anti-clusters
export CountAllClusters=2       ## = 0: only physical clusters are counted (unphysical clusters are counted as single baronys), 
                                ## = 1: all clusters A > 7 are counted as clusters
				## = 2: both options 0 and 1 are executed one after the other

## Options for freeze-out coordinates
export PhqmdWithFreeze=1        ## = 0: PHQMD code without freeze-out coordinates, = 1: PHQMD code with freeze-out coordinates
# For option PhqmdWithFreeze = 1:
export WriteUnigen=0            ## = 1: UniGen-file (including freeze-out position) is written
export WriteEventFreeze=1       ## = 1: EventFreeze-file (including freeze-out posiition and momentum) is written

## Delete & Merge
export DeleteFiles=0            ## = 1: deletes all output files, but the detector-input [dataset].phqmd.root
export RunHadd=0                ## = 1: add 100 jobs into one file
export RunHaddSim=0             ## = 1: add 10 jobs in one file for detector simulations

###############################################################
## Parameters for sbatch ######################################

export firstJob=1
export lastJob=1000
export array=$firstJob-$lastJob
export last_hadd=$(printf "%0.0f" "$(echo " $lastJob / 100  " | bc -l)")
export array_hadd=1-$last_hadd
export last_hadd_simcbm=$(printf "%0.0f" "$(echo " $lastJob / 10  " | bc -l)")
export array_hadd_simcbm=$firstJob-$last_hadd_simcbm
export time="5-23:59:59"
export partition="long"
export ram="8G"

###############################################################
## Parameters for inputPHSD (more parameters availabe, but should usually not be changed)

export SYSTEM=auau              ## specification below for auau, pbpb, aupb, aupt - other systems need to be added below
export ENERGY=4.93              ## center of mass energy - !for sqrt(s) > 7 GeV setting for kinetic deuterons in PHQMD code should be changed! 
export ELAB=$(echo " $ENERGY * $ENERGY / (2 * 0.938) - (2 * 0.938)" | bc -l)  ##  Lab energy per nucleon (needed for PHQMD)
export NUM=100                  ## number of parallel events (NUM=100 for 3 GeV & 4.9 GeV required)
export ISUBS=10                 ## number of subsequent runs
export nEvents=$(echo " $NUM * $ISUBS " | bc -l)

export IMPACTPARAMETER_MIN=0.0  ## minimal impact parameter in fm
export IMPACTPARAMETER_MAX=15.0 ## maximum impact parameter in fm ! Max impact parameter must not be > 15.0 fm !

export IGLUE=1                  ## =1 with partonic QGP phase (PHSD mode); =0 - HSD mode ! Needs to be set to = 0 for 3 GeV
export EOS=0                    ## EoS; =0: hard EOS without M.D.I; =1 soft EoS; =2 soft EoS with mom. dependence 
export ICLUSTER=1               ## enable or disable CLUSTER output

##############################################################
#############  Directories and files #########################
    
export version_phqmd=PHSD-PHQMD
export LOCATION=/lustre
export DIR=$LOCATION/cbm/users/$USER
export PHQMDDIR=$DIR/$version_phqmd
export ROOT_SOURCE=/cvmfs/fairsoft.gsi.de/debian11/fairsoft/nov22p1/bin
export UNIGEN_SOURCE=$DIR/unigen
#export CBMROOT_DIR=$DIR/cbmroot/build
export USE_CBMROOT=0           ## =1: if cbmroot should be sourced instead of root and UniGen
export OUTDIR=$DIR/mc/$version_phqmd/$SYSTEM"_"$ENERGY"GeV_"$IMPACTPARAMETER_MAX"fm_Num"$NUM"xSub"$ISUBS"_EoS"$EOS

##############################################################
############### Specify collision system #####################

if [[ $SYSTEM == "auau" ]] ; then
    export MASSTA=197
    export MSTAPR=79
    export MASSPR=$MASSTA
    export MSPRPR=$MSTAPR
fi

if [[ $SYSTEM == "pbpb" ]] ; then
    export MASSTA=208
    export MSTAPR=82
    export MASSPR=$MASSTA
    export MSPRPR=$MSTAPR
fi

if [[ $SYSTEM == "aupb" ]] ; then
    export MASSTA=208
    export MSTAPR=82
    export MASSPR=197
    export MSPRPR=79
fi

if [[ $SYSTEM == "aupt" ]] ; then                                                                                                           
    export MASSTA=195                                                                                                                       
    export MSTAPR=78                                                                                                                        
    export MASSPR=197                                                                                                                       
    export MSPRPR=79                                                                                                                        
fi




