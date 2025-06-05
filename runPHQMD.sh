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

if [ "$CountAllClusters" == 1 ] || [ "$CountAllClusters" == 2 ]; then
    export FOLDER_CL=allclusters
else
    export FOLDER_CL=smallclusters
fi

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

if [ "$RunStabilisation" == 1 ]; then
    jobname=stab
    echo "Run Stabilisation"
    cp $SCRIPTDIR/791to891.f $OUTDIR
    export script_stab=791to891.exe
    cp $SCRIPTDIR/$script_stab $OUTDIR
fi

if [ "$MakeDetectorInput" == 1 ]; then
    jobname=conv
    echo "Convert PHQMD output into detector input"
    if [ "$PhqmdWithFreeze" == 0 ]; then
	    export script_convert="convert_phqmd_detector_unigen.C"
    fi
    if [ "$PhqmdWithFreeze" == 1 ]; then
	    export script_convert="convert_phqmd_detector_unigen_freezeout.C"
    fi
    mkdir -p $OUTUNIGEN
    mkdir -p $OUTUNIGEN/$FOLDER_CL
    echo "UniGen files will be written to: ${OUTUNIGEN}/${FOLDER_CL}"  
    cp $SCRIPTDIR/$script_convert $OUTUNIGEN/$FOLDER_CL
    cp $SCRIPTDIR/cluster_table.dat $OUTUNIGEN/$FOLDER_CL
    if [ "$USE_CBMROOT" == 0 ]; then 
	 cp $UNIGEN_SOURCE/rootlogon.C $OUTUNIGEN/$FOLDER_CL
    fi
    if [ "$CountAllClusters" == 2 ]; then
	 export FOLDER_CL=smallclusters
	 mkdir -p $OUTUNIGEN/$FOLDER_CL                             
         echo "UniGen files for small clusters will be written to: ${OUTUNIGEN}/${FOLDER_CL}" 	 
	 cp $SCRIPTDIR/$script_convert $OUTUNIGEN/$FOLDER_CL                                                                                  
	 cp $SCRIPTDIR/cluster_table.dat $OUTUNIGEN/$FOLDER_CL                                                                                
         if [ "$USE_CBMROOT" == 0 ]; then
      	      cp $UNIGEN_SOURCE/rootlogon.C $OUTUNIGEN/$FOLDER_CL
	 fi	
    fi
fi

if [ "$DeleteFiles" == 1 ]; then
    echo "Delete files, only keep detector input"
fi 

if [ "$RunHadd" == 1 ]; then
	jobname=hadd
	echo "Run hadd" 
	export INDIR=$OUTUNIGEN/$FOLDER_CL
	export OUTDIR=$INDIR/hadd
	echo "Output will be written to: ${OUTDIR}"
	mkdir -p $OUTDIR
	export LOGDIR=$OUTDIR"/log"
	mkdir -p $LOGDIR
        if [ "$ConvertAntiClusters" == 1 ]; then
		export file=phqmd.root
	else
		export file=phqmd_noanti.root
	fi
	export array=$array_hadd
	export batchfile=$SCRIPTDIR/batch_hadd_file.sh
fi

if [ "$RunHaddSim" == 1 ]; then
	jobname=haddsim
	echo "Run hadd for simcbm"
      	export INDIR=$OUTUNIGEN/$FOLDER_CL
	export OUTDIR=$INDIR/haddsimcbm
	echo "Output will be written to: ${OUTDIR}"
	mkdir -p $OUTDIR
	export LOGDIR=$OUTDIR"/log"
	mkdir -p $LOGDIR
	if [ "$ConvertAntiClusters" == 1 ]; then                                                                                                                                  
		export file=phqmd.root                                
	else                                                          
		export file=phqmd_noanti.root    
	fi                                                                                                                                                                        
	export array=$array_hadd_simcbm                                                                                                                   
	export batchfile=$SCRIPTDIR/batch_hadd_file.sh                                                                                                                     
fi 

sbatch --job-name=${jobname} --partition=${partition} --time=${time} --array=${array} -D $LOGDIR -o %a_%A.out.log -e %a_%A.err.log --export=ALL -- $batchfile

if [ "$CountAllClusters" == 2 ]; then
   if [ "$RunHadd" == 1 ] || [ "$RunHaddSim" == 1 ]; then
	  export FOLDER_CL=smallclusters
          export INDIR=$OUTUNIGEN/$FOLDER_CL 
          
	  if [ "$RunHadd" == 1 ]; then 
	       jobname=hadd
	       export OUTDIR=$INDIR/hadd
               export array=$array_hadd
	  fi
          if [ "$RunHaddSim" == 1 ]; then 
	       jobname=haddsim
               export OUTDIR=$INDIR/hadd
               export array=$array_hadd_simcbm
	  fi 

          echo "Output for small clusters will be written to: ${OUTDIR}"                                                                                      
	  mkdir -p $OUTDIR                                                                                                                 
	  export LOGDIR=$OUTDIR"/log"                                                                                                      
	  mkdir -p $LOGDIR 
	  export batchfile=$SCRIPTDIR/batch_hadd_file.sh	  

	  sbatch --job-name=${jobname} --partition=${partition} --time=${time} --array=${array} -D $LOGDIR -o %a_%A.out.log -e %a_%A.err.log --export=ALL -- $batchfile
     fi
fi
