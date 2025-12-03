#!/bin/sh
export runscript="runinfo.sh"
export SCRIPTDIR=$(pwd)
export PATH=$PATH:$SCRIPTDIR
export batchfile="batch_runPHQMD.sh"

. $SCRIPTDIR/$runscript

echo "Lab energy is set to: ${ELAB} AGeV"

mkdir -p $OUTDIR
export LOGDIR=$OUTDIR"/log"
mkdir -p $LOGDIR

export OUTUNIGEN=$OUTDIR/unigen

cp $SCRIPTDIR/$runscript $OUTDIR
cp $SCRIPTDIR"/runPHQMD.sh" $OUTDIR
cp $SCRIPTDIR/$batchfile $OUTDIR

if [ "$RunPhqmdCode" == 1 ]; then
    jobname=phqmd    
    echo "Run PHQMD code"
    echo "PHQMD output will be written to: ${OUTDIR}"
    cp $PHQMDDIR/phqmd_info $OUTDIR
fi

cd $SCRIPTDIR

if [ "$MakeDetectorInput" == 1 ]; then
    jobname=conv
    echo "Convert PHQMD output into detector input"
    if [ "$PhqmdWithFreeze" == 0 ]; then
	    export script_convert="convert_phsd_detector_unigen.C"
    fi
    if [ "$PhqmdWithFreeze" == 1 ]; then
	    export script_convert="convert_phsd_detector_unigen_freezeout.C"
    fi
    mkdir -p $OUTUNIGEN

    echo "UniGen files will be written to: ${OUTUNIGEN}"  
    cp $SCRIPTDIR/$script_convert $OUTUNIGEN
    cp $SCRIPTDIR/cluster_table.dat $OUTUNIGEN
    if [ "$USE_CBMROOT" == 0 ]; then 
	 cp $UNIGEN_SOURCE/rootlogon.C $OUTUNIGEN
    fi
fi

if [ "$DeleteFiles" == 1 ]; then
    echo "Delete files, only keep detector input"
fi 

if [ "$RunHadd" == 1 ]; then
	jobname=hadd
	echo "Run hadd" 
	export INDIR=$OUTUNIGEN
	export OUTDIR=$INDIR/hadd
	echo "Output will be written to: ${OUTDIR}"
	mkdir -p $OUTDIR
	export LOGDIR=$OUTDIR"/log"
	mkdir -p $LOGDIR
	export file=phsd.root
	export array=$array_hadd
	export batchfile=$SCRIPTDIR/batch_hadd_file.sh
fi

if [ "$RunHaddSim" == 1 ]; then
	jobname=haddsim
	echo "Run hadd for simcbm"
      	export INDIR=$OUTUNIGEN
	export OUTDIR=$INDIR/haddsimcbm
	echo "Output will be written to: ${OUTDIR}"
	mkdir -p $OUTDIR
	export LOGDIR=$OUTDIR"/log"
	mkdir -p $LOGDIR                                                                                                                            
	export file=phsd.root                                                                                                                                                                                                  
	export array=$array_hadd_simcbm                                                                                                                   
	export batchfile=$SCRIPTDIR/batch_hadd_file.sh                                                                                                                     
fi 

sbatch --job-name=${jobname} --partition=${partition} --time=${time} --array=${array} -D $LOGDIR -o %a_%A.out.log -e %a_%A.err.log --export=ALL -- $batchfile
