#!/bin/bash

. $ROOT_SOURCE/thisroot.sh
#. $CBMROOT_DIR/config.sh -a

if [ "$RunHadd" == 1 ]; then
	XXX=$(printf "%03d" "$SLURM_ARRAY_TASK_ID")
	XINDEX=$(echo " ($SLURM_ARRAY_TASK_ID - 1) * 100 " | bc -l)
	lastJob=101
	outfile=$OUTDIR/$XXX.hadd1.$file 
fi

if [ "$RunHaddSim" == 1 ]; then
	XXX=$(printf "%04d" "$SLURM_ARRAY_TASK_ID")
	XINDEX=$(echo " ($SLURM_ARRAY_TASK_ID - 1) * 10 " | bc -l)
	lastJob=11
	outfile=$OUTDIR/$XXX.$file 	
fi	

echo "XXX=${XXX}"

iJob=1
allinFiles=""

while [ "$iJob" != "$lastJob" ]; do
    number=$(echo " $XINDEX + $iJob " | bc -l)
    IIIII=$(printf "%05d" "$number")

    inFile=$INDIR/$IIIII.$file
    
    echo "inFile: ${inFile}"
    if [ -e $inFile ]; then
	echo "... is there"
	allinFiles="${allinFiles} ${inFile}"
    else
	echo "... is missing"
    fi
    
    let iJob=$iJob+1;
done

cd $INDIR

rm $outfile
root -l -b -q "rootlogon.C"
hadd -k ${outfile} ${allinFiles}

cd $SCRIPTDIR
