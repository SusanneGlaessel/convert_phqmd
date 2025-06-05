# RUNNING & CONVERTING PHQMD

### 3 steps are required:

1) run the PHQMD code
2) stabilize the cluster output
3) convert the PHQMD-output into detector-input (UniGen format)
optional:
4) delete all files, keep only the detector input
5) merge output files - 100 or 10 outputfiles (for detector simulations) into 1 file
All 4 steps are performed with the batch-script. Single steps can be switched on or off. 

### The following files are needed (besides the PHQMD code):

- runinfo.sh
- runPHQMD.sh
- batch_runPHQMD.sh
- 791to891.exe*
- convert_phqmd_detector_unigen.C / convert_phqmd_detector_unigen_freeze.C***** 
- cluster_table.dat***

### Only runinfo.sh needs to be modified:

- selection of steps
- option for cluster conversion mode**
- choose if only clusters OR clusters & anti-clusters are converted****
- options for freeze-out coordinates*****
- select input information for PHQMD, e.g system, energy, number of parallel events,
impactparamter etc. 
- location information

### Run all 3 steps:

	. runPHQMD.sh

### Comments:

\* Stabilisation: 791to891.exe converts fort.791 (fort.781) into fort.891 (fort.881) which
considers a clusters as stable when the cluster-baryons are freezed out (and if it has a 
negative binding energy).

** Conversion: convert_phqmd_detector_unigen.C creates the detector input from the
PHQMD output. It converts the PHQMD output-files after stabilization to the UniGen
format. 
In PHQMD clusters / anticlusters are recognised independently of their physical
existence. The cluster-baryons are listed separately in the outputfile. This routine
builds clusters from the single baryons based on their cluster-ID and identifies physical
clusters according to the cluster_table.root. Baryons from unphysical clusters (eg. p-p)
are counted as single baryons. With option "CountAllClusters" all clusters with size A > 7 
are counted as clusters independent of their physical existence.

*** The cluster_table.dat contains the information about physical clusters, their baryon 
content and branching ratio. It is required to perform the conversion and can be easily 
modified/extended. Each line contains the following information of the respective cluster:

	name / pdgcode / number of protons / number of neutral baryons / number of Lambdas 
	number of Simga0 / branching ratio

**** PHQMD writes baryons and anti-baryons into two separate files. The conversion
of anti-baryons is optional and can be switched off.

***** PHQMD output files can be written with freeze-out coordinates. This needs to be changed in 
the PHQMD-code itself. In order to read the ouput with freeze-out coordinates PhqmdWithFreeze 
needs to be set to 1. The ouput-format can either be UniGen (which does not include freeze-out 
mometum) are EventFreeze(which includes the freeze-out mometum).

For questions, please contact glaessel@ikf.uni-frankfurt.de

