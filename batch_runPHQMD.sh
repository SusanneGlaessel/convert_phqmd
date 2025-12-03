#!/bin/bash

INDEX=$SLURM_ARRAY_TASK_ID
XXXXX=$(printf "%05d" "$INDEX")
echo "SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}"
echo "Output directory is: ${OUTDIR}"

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
 "$MASSPR",       MASSPR: projectile mass (=0 for pi- beam)
 "$MSPRPR",        MSPRPR: protons in projectile
 "$ELAB",      ELAB:   kinetic energy per nucleon in Lab frame ! Tkin=21300 (srt=200 GeV) RHIC, =13433049 (5 TeV) = 26120000 (7 TeV) 
 "$IMPACTPARAMETER_MIN",       BMIN:   minimal impact parameter in fm
 "$IMPACTPARAMETER_MAX",       BMAX:   maximal impact parameter in fm
 0.5,       DBIMP: impact parameter step in fm (used only if IBweight_MC=0)
 "$NUM",       NUM:  number of parallel ensambles (events) (>10 - see manual!)
 "$ISUBS",         ISUBS:  number of subsequent runs
 "$ISEED",    ISEED:  ANY uneven INTEGER number
 "$IGLUE",         IGLUE: =1 with partonic QGP phase (PHSD mode); =0 - HSD mode 
 145.0,     FINALT: final time of calculation in fm/c (if set =0, it will be computed in the PHSD)
 10,        ILOW: output level (default=10)
 0,         Idilept: =0 no dileptons; =1 electron pair; =2  muon pair
 0,         ICQ: =0 free rho's, =1 dropping mass, =2 broadening, =3 drop.+broad.
 0,         IHARD: =1 with charm and bottom; =0 - without 
 1,         iDQPM: =1 with DQPM(T,muB) for QGP, =0 as in PHSD4.X - old DQPM(T) (default=1)
 1,         IBweight_MC: =0 constant step in B =DBIMP; =1 choose B by Monte-Carlo ISUBS times in [Bmin,Bmax]
 0,         IUSER: =1 for general users : use default /optimized settings; = 0 for PHSD team
 1,         INUCLEI: =1 reactions with deuterons=1 including kinetic reactions with deuterons  
 1,         IPHQMD=1: propagation with QMD dynamics; =0 with HSD/PHSD dynamics !-->>> flags below for PHQMD ONLY:
 0,         ISACA: enable or disable SACA output (note: SACA or MST is controled by iflagsaca)
 0.0,       TSACA: starting time for SACA
 0.0,       DTSACA, time step for SACA calculations
 0,         NTSACA, Number of SACA timesteps ( =>  tmax=tsaca+dtsaca*ntsaca must be < FINALT !)
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
 0,         IOUT:  default=0                                                                                           
 0,          IRELAX: default=0  !<<<---- end QMD flags ---------------
 0,       IANTIPROT: =-1 for antiproton beam (projectile), =0 for all other projectiles (default=0)
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

    phsdFile="phsd.dat"
    cd $OUTDIR/$XXXXX

    if ! [ -e $phsdFile ]; then                         
	    echo "Error: ${phsdFile} is missing"
	    exit 1
    fi

    if [ "$PhqmdWithFreeze" == 0 ]; then	
	 cd $OUTUNIGEN
	 root -l -b -q  "$script_convert(\"$OUTDIR\",\"$XXXXX\",$firstevent)"	      
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
       
       	cd $OUTUNIGEN
	root -l -b -q  "$script_convert(\"$OUTDIR\",\"$XXXXX\",$firstevent,$IsWriteUnigen,$IsWriteEventFreeze)"
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
