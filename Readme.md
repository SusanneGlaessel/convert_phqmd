### CONVERT PHQMD

## General information
convert_phqmd converts the PHQMD output files (after stabilisation)
into the UniGen format.

Branch **PHQMD52 _Winn_clean** is compatible with PHQMD version
PHQMD52-Winn_2023. It does not include the scripts to  run the PHQMD-code & the stabilisation .

For the scripts to run the PHQMD-code & the stabilisation in addition, please use  **PHQMD52 _Winn**.

## Requirements

### Root
https://root.cern/install/

Root with c++17 standard is recommended.

### Unigen
https://github.com/FairRootGroup/UniGen

## Installation

Clone

	git clone --branch PHQMD52_Winn_clean git@github.com:SusanneGlaessel/convert_phqmd

Set rootsource

	source path_to_root/bin/thisroot.sh

Load unigen libraries

	cd convert_phqmd
	root -l rootlogon.C
	
The *path\_to\_installation* in rootlogon.C needs to be replaced by the acutal
location of UniGen:

	LINUX: gSystem->Load("path_to_unigen_installation/install/lib/l libunigen.so");
	MAC:
	gSystem->Load("path_to_unigen_installation/install/lib/libunigen.dylib");

## Inputfiles

The following output files from PHQMD are needed:
- inputPHSD
- phsd.dat
- fort.891
- fort.881 (optional)

## Clustertable

In PHQMD clusters / anticlusters are recognised independently of their physical
existence. This routine identifies physical clusters according to the cluster_table.dat. Baryons from unphysical clusters (eg. p-p)
are counted as single baryons. 

The cluster_table.dat can be easily modified/extended. Each line contains the following
information of the respective cluster:
 
	name / pdgcode / number of protons / number of neutral baryons /
	number of Lambdas / number of Simga0 / branching ratio
	
## Run

Run macro:
 
	root -l convert_phqmd.C

 Input directory and options can be managed by the settings:

 - ConvertAllClusters:
	With option ConvertAllClusters = kTRUE all clusters with size A > 7 
are counted as clusters independent of their physical existence.

- ConvertAnti:
   PHQMD writes baryons and anti-baryons into two separate files. The conversion
of anti-baryons is optional and can be switched off with ConvertAnti =
kFALSE.

 - FreezeCoords:
    PHQMD output files can be written with freeze-out coordinates. This needs to be changed in 
the PHQMD-code itself. In order to read the ouput with freeze-out
coordinates FreezeCoords = kTRUE needs to be set.
     The freeze-out momtum is not written into UniGen. To create a
    format including the freeze-out momentum, use branch **PHQMD52 _Winn**.

***** 

For questions, please contact glaessel@ikf.uni-frankfurt.de

