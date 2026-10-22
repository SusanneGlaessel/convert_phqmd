# RUNNING & CONVERTING PHQMD without clusters

## 2 steps need to be performed:

1) run the PHQMD code
2) convert the PHQMD-output into detector-input (UniGen format)

optional:
3) delete all files, keep only the detector input
4) merge output files - 100 into 1 file
All steps are performed with the batch-script. Single steps can be
switched on or off.

## Requirements (besides the PHQMD code)

### Root
https://root.cern/install/

Root with c++17 standard is recommended.

### UniGen
https://github.com/FairRootGroup/UniGen

The *path\_to\_unigen\_installation* in macro/rootlogon.C needs to be replaced by the acutal
location of UniGen:

	LINUX: gSystem->Load("path_to_unigen_installation/install/lib/libunigen.so");
	MAC:   gSystem->Load("path_to_unigen_installation/install/lib/libunigen.dylib");

## The following files are needed for running and converting:

- runinfo.sh
- runPHQMD.sh
- batch_runPHQMD.sh
- convert_phsd_detector_unigen.C / convert_phsd_detector_unigen_freeze.C***

## Only runinfo.sh needs to be modified:

- selection of steps
- options for freeze-out coordinates***
- select input information for PHQMD, e.g system, energy, number of parallel events,
impactparamter etc. 
- location information

## Run all 2 steps:

	. runPHQMD.sh

## Comments:

** Conversion: convert_phqmd_detector_unigen.C creates the detector input from the
PHQMD output. It converts the PHQMD output-file to the UniGen
format. 

*** PHQMD output files can be written with freeze-out coordinates. This needs to be changed in 
the PHQMD-code itself. In order to read the ouput with freeze-out coordinates PhqmdWithFreeze 
needs to be set to 1. The ouput-format can either be UniGen (which does not include freeze-out 
mometum) are EventFreeze(which includes the freeze-out mometum).

For questions, please contact glaessel@ikf.uni-frankfurt.de

