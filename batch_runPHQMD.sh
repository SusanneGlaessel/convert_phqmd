#!/bin/bash

INDEX=$SLURM_ARRAY_TASK_ID
XXXXX=$(printf "%05d" "$INDEX")
echo "SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}"
echo "Output directory is: ${OUTDIR}"

## Time settings for PHQMD & Stabilisation
TSACA=10.0               ## TSACA: starting time for SACA 
NTSACA=26                ## NTSACA, Number of SACA timesteps ( =>  tmax=tsaca+dtsaca*ntsaca must be < FINALT !)
DTSACA=5.0               ## DTSACA, time step for SACA calculations (default=25.)
FINALT=145.0             ## final time of calculation in fm/c

##################################################################################
################### Run PHQMD ####################################################
##################################################################################

if [ "$RunPhqmdCode" == 1 ]; then
    
    echo "Run PHQMD"
    echo "PHQMDDIR = ${PHQMDDIR}"

    cd $OUTDIR
    mkdir -p $XXXXX
    cd $XXXXX

    ## Create ISEED
    year=`date +%Y`
    day=`date +%d`
    month=`date +%m`
    hour=`date +%H`
    min=`date +%M`
    sec=`date +%S`
    ISEED=$(printf "%.f" "$(echo "($sec * $year + $SLURM_ARRAY_TASK_ID * $day * $hour * $min + $month) *2 + 1" | bc -l)")
    echo ""$ISEED"" > "ISEED.txt"
    echo "ISEED=${ISEED}"

    ## Copy PHQMD-code into outputfolder
    ln -s $PHQMDDIR/HSDINPUT .
    ln -s $PHQMDDIR/CC_INPUT .
    ln -s $PHQMDDIR/phqmd .
    ln -s $PHQMDDIR/iso.data .
    ln -s $PHQMDDIR/mass.inp .

    ## Create inputPHSD
    echo " "$MASSTA",       MASSTA: target mass / au=197 / pb=208
 "$MSTAPR",        MSTAPR: protons in target / au=79  / pb=82
 "$MASSPR",       MASSPR: projectile mass
 "$MSPRPR",        MSPRPR: protons in projectile
 "$ELAB",      ELAB:   4060000. Lab energy per nucleon (e.g.: 200GeV, b=2.2fm: 5%, t(CPU)~4h)
 "$IMPACTPARAMETER_MIN",       BMIN:   minimal impact parameter in fm
 "$IMPACTPARAMETER_MAX",       BMAX:   maximal impact parameter in fm
 0.5,       DBIMP: impact parameter step in fm (used only if IBweight_MC=0)
 "$NUM",       NUM:    number of parallel events(>10)
 "$ISUBS",         ISUBS:  number of subsequent runs
 "$ISEED",    ISEED:  ANY uneven INTEGER number
 "$IGLUE",         IGLUE: =1 with partonic QGP phase (PHSD mode); =0 - HSD mode 
 "$FINALT",     FINALT: final time of calculation in fm/c
 10,        ILOW: output level (default=10)
 0,         Idilept: =0 no dileptons; =1 electron pair; =2  muon pair
 0,         ICQ: =0 free rho's, =1 dropping mass, =2 broadening, =3 drop.+broad.
 0,         IHARD: =1 with charm and bottom; =0 - without 
 1,         iDQPM: =0 PHSD4.X, =1 PHSD5.X with (T,muB) for QGP
 1,         IBweight_MC: =0 constant step in B =DBIMP; =1 choose B by Monte-Carlo ISUBS times in [Bmin,Bmax]
 0,         IUSER: =1 for general users : use default /optimized settings; = 0 for PHSD team
 1,         INUCLEI  =1 reactions with deuterons 
 1,         IPHQMD=1: propagation with QMD dynamics; =0 with HSD/PHSD dynamics
 "$ICLUSTER",         ISACA: enable or disable SACA output (note: SACA or MST is controled by iflagsaca)
 "$TSACA",      TSACA: starting time for SACA
 "$DTSACA",       DTSACA, time step for SACA calculations
 "$NTSACA",        NTSACA, Number of SACA timesteps ( =>  tmax=tsaca+dtsaca*ntsaca must be < FINALT !)
 0,         iflagsaca: =1 SACA(=FRIGA) analysis(+MST), =0 MST, no SACA --!!!! FRIGA
 0,         tageyuk=1  ! if Yukawa potential requested in SACA       
 0,         tageasy=1  ! if asymmetry energy requested in SACA 
 0,         tagepair=1 ! tructure effects (pairing,...) requested in SACA 
 0,         tagecoul=0 ! activates the coulomb energy for fragment selection
 23.3,      vasy0=23.3 ! asymmetry potential energy at normal density (MeV) in SACA 
 0.0,       eta_pairing=0.0 ! pairing potential exponant (0.->only forbids unbound isotopes, 1., 0.65, 0.35 or 0.25) in SACA 
 "$EOS",         iqmdeos  ! EoS for QMD option IPHQMD=1 ; =0: hard EOS without M.D.I; =1 soft EoS; =2 soft EoS with mom. dependence
 3,         IFLAG_Res_SACA ! =1 include ALL resonances with their decay to SACA; =2 - only nucleons; =3 nucleons and hyperons
 0,         IfragWigDen  ! =0: no; =1 yes, light clusters formation according to the Wigner density
 3,         ILOC =3 location of HSD data files                                                                
 0,         IOUT =0                                                                                           
 0,         IRES
 " > "inputPHSD"

    time ./phqmd

    mv $OUTDIR/$XXXXX/inputPHSD $OUTDIR/

    ## Remove PHQMD-Code
    rm -r $OUTDIR/$XXXXX/HSDINPUT
    rm -r $OUTDIR/$XXXXX/CC_INPUT
    rm $OUTDIR/$XXXXX/phqmd
    rm $OUTDIR/$XXXXX/iso.data
    rm $OUTDIR/$XXXXX/mass.inp

    ## Remove not needed outputfiles
    rm -r $OUTDIR/$XXXXX/OUTPUT
    rm $OUTDIR/$XXXXX/testDQPM.txt
    rm $OUTDIR/$XXXXX/fort.888
    rm $OUTDIR/$XXXXX/fort.80
    rm $OUTDIR/$XXXXX/fort.790
    rm $OUTDIR/$XXXXX/fort.780
    
    cd $SCRIPTDIR

fi

##################################################################################
################### Stabilisation routine for fort.791 & fort.781 ################
##################################################################################

if [ "$RunStabilisation" == 1 ]; then
    
    phsdFile="phsd.dat"
    cd $OUTDIR/$XXXXX

    if ! [ -e $phsdFile ]; then                         
	    echo "Error: ${inputFile} is missing"
	    exit 1
    fi

    echo "Run stabilisation routine for baryons from fort.791"

    ## Create parameter file
    echo ""$NUM",          number of parallel ensembles
"$NTSACA",           number of timesteps
"$TSACA",         initial time to start SACA output
"$DTSACA",          size of timestep
fort.791,     input filename
fort.891,     output filename
" > "parameters.txt"

    cd $OUTDIR/$XXXXX
    cp $SCRIPTDIR/$script_stab .
    ./$script_stab
    rm parameters.txt

    if [ "$ConvertAntiClusters" == 1 ]; then
	
	echo "Run stabilisation routine for anti-baryons from fort.891"

	## Create parameter file
	echo ""$NUM",          number of parallel ensembles
"$NTSACA",           number of timesteps
"$TSACA",         initial time to start SACA output
"$DTSACA",          size of timestep
fort.781,     input filename
fort.881,     output filename
" > "parameters.txt"

	cd $OUTDIR/$XXXXX
	./$script_stab
	rm $script_stab
	rm parameters.txt
    fi
fi

##################################################################################
######## Conversion of PHQMD output into detector input (unigen format) ##########
##################################################################################

if [ "$MakeDetectorInput" == 1 ]; then
    
    echo "Start conversion of PHQMD output files into unigen format"

    export firstevent=$(echo "($SLURM_ARRAY_TASK_ID-1)*$ISUBS*$NUM" | bc -l)
    echo "firstevent is: ${firstevent}"

    if [ "$USE_CBMROOT" == 1 ]; then
         . $CBMROOT_DIR/config.sh -a
    else
	 . $ROOT_SOURCE/thisroot.sh
    fi

    inputFile="fort.891"
    cd $OUTDIR/$XXXXX

    if ! [ -e $inputFile ]; then
	echo "Error: ${inputFile} is missing"
	exit 1
    fi
    
    if [ "$ConvertAntiClusters" == 1 ]; then
	
	inputFile="fort.881"
	cd $OUTDIR/$XXXXX

	if ! [ -e $inputFile ]; then
	    echo "Error: ${inputFile} is missing"
	    exit 1
	fi
	
	export IsConvertAntiClusters=kTRUE
     
     else
        export IsConvertAntiClusters=kFALSE
    fi

    if [ "$CountAllClusters" == 1 ] || [ "$CountAllClusters" == 2 ] ; then
	export IsCountAllClusters=kTRUE  
        export FOLDER_CL=allclusters	
    else                                                                   
	export IsCountAllClusters=kFALSE
	export FOLDER_CL=smallclusters
    fi

    if [ "$PhqmdWithFreeze" == 0 ]; then	
	 cd $OUTUNIGEN/$FOLDER_CL
	 root -l -b -q  "$script_convert(\"$OUTDIR\",\"$XXXXX\",$firstevent,$IsCountAllClusters,$IsConvertAntiClusters)"
         
	 if [ "$CountAllClusters" == 2 ] ; then
              export IsCountAllClusters=kFALSE
              export FOLDER_CL=smallclusters
	      cd $OUTUNIGEN/$FOLDER_CL  
              root -l -b -q  "$script_convert(\"$OUTDIR\",\"$XXXXX\",$firstevent,$IsCountAllClusters,$IsConvertAntiClusters)"
	 fi	      
    fi	
    
    if [ "$PhqmdWithFreeze" == 1 ]; then
	if [ "$WriteUnigen" == 1 ]; then
            export IsWriteUnigen=kTRUE
        else
            export IsWriteUnigen=kFALSE
	fi

	if [ "$WriteEventFreeze" = 1 ]; then
	    export IsWriteEventFreeze=kTRUE
        else
	    export IsWriteEventFreeze=kFALSE
	fi
       
       	cd $OUTUNIGEN/$FOLDER_CL
	root -l -b -q  "$script_convert(\"$OUTDIR\",\"$XXXXX\",$firstevent,$IsCountAllClusters,$IsConvertAntiClusters,$IsWriteUnigen,$IsWriteEventFreeze)"
        
	if [ "$CountAllClusters" == 2 ] ; then
             export IsCountAllClusters=kFALSE
	     export FOLDER_CL=smallclusters
	     cd $OUTUNIGEN/$FOLDER_CL
	     root -l -b -q  "$script_convert(\"$OUTDIR\",\"$XXXXX\",$firstevent,$IsCountAllClusters,$IsConvertAntiClusters,$IsWriteUnigen,$IsWriteEventFreeze)"
	fi
    fi	
fi

##################################################################################
######## Delete all files, but the detector input ##########
##################################################################################

if [ "$DeleteFiles" == 1 ]; then
    echo "Delete files, only keep detector input"
    rm -r $OUTDIR/$XXXXX
fi 
    
cd $SCRIPTDIR
