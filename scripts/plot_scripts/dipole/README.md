# Anisotropy_Apj_2022


Code to reproduce the results in the paper "Arrival Directions of Cosmic Rays above 32 EeV from Phase One of the Pierre Auger Observatory".


##Required externals

This code requires some external software to work:

-[] [Cern ROOT v6.xx] {https://root.cern} is needed in all of the analyses to create and fill histograms, functions, and use formulas. 

-[] [HEALPix v3.xx] {https://healpix.sourceforge.io/downloads.php} is a sky pixelization software needed in the blind and catalog-based searches. 

-[] Python3 or Phython 2.7 is necessary for recursion scripts and figure plotting.



##Code Structure

The code itself is divided into four main parts.

-[] The common utilities and methods in the Utilities folder.

-[] The data and exposure calculation routine in the Data folder.

-[] The main analyses in the Targeted_and_Blind and Catalog_Based folders.

-[] The figure plotting scripts in the Visuals folder.


The order of execution of the code is essential: exposure calculation in Data >> Core analyses >> Figure plotting in Visuals.


#Utilities

The Utilities folder contains methods that are of common use throughout all the analyses, like angular distance calculations, tools for healpix maps and cosmetics. It is automatically compiled by any of the Makefiles.


#Data

The Data folder contains only one executable exposure.exe that needs no additional arguments and produces two output files, exposure.fits and exposure.root needed by the analyses.
The DataPath.h file contains the declaration of the dataset file to use in all of the analyses. The time_exposure.dat file provides the evolution of exposure (in km^2 sr yr) with time (in human readable format).

Running ./exposure.exe is REQUIRED for running the other analyses. Typical running time: minutes.


#Targeted_and_Blind

This folder contains three subfolders for the Blind, Autocorrelation and Targeted (Galactic plane, Galactic center, Supergalactic plane, Centaurus region) analyses and the dedicated utilities.
All folders have the same structure, with two executables, 1_pvalue.exe for the local p-value calculation, and 2_penalization.exe which outputs the post trial p-value.

-[] In the Blind subfolder, the 1_pvalue.exe executable requires the number of isotropic events to be simulated as argument; 2_penalization.exe requires the number of isotropic skies and the number of isotropic events. Examples: ./1_pvalue.exe 1000000; ./2_penalization.exe 10000 1000000. Typical running time: 10 minutes local p-value. 1 hour every 100 simulations penalization.

-[] In the Autocorrelation subfolder, both executables require as argument the number of isotropic skies to be simulated. The results of the local p-value are stored in ROOT files. Examples: ./1_pvalue.exe 100000; ./2_penalization.exe 100000. Typical running time: 1 hour every 1000-1500 simulations.

-[] In the Targeted subfolder, both executables require as arguments the type of analysis ("gc", "gp", "sgp", "ca") and the number of isotropic skies to simulate. The 2_penalization.exe exectuable requires also an additional argument that can be 0 if one wants to penalize for the scan in angle only, and any other integer if one wants to perform a full penalization. Examples: ./1_pvalue.exe gc 1000000; ./2_penalization.exe ca 1000000 1. Typical running time: 1 hour every 10000-15000 simulations.

All the executables, if ran without specifying the number of simulations, default to a base value specified in the utils.


#Catalog_Based

This folder contains the code for the catalog based analysis, its dedicated utilities and data.

The subfolder Catalogs contains a C/ROOT executable, propagate.exe, a python script, propagate.py, a ROOT file, AnaCRP3.root, and two subfolders, Multiwavelength, that contains the input of the executable, and ModelsUHECR that contains the output. Running propagate.py populates the output subfolder. Executing this python script is sufficient but REQUIRED to run the analysis. Typical running time: 1 hour.

The executable 1_pre_trial.exe stores results in ROOT files; 2_penalization.exe produces post trial p-values as terminal output. Typical running time: 30 minutes.


#Visuals

This folder contains two Python scripts and a C/ROOT executable. skymaps.py produces figures 1, 8 and 10 in the paper; show_figures.py produces figures 2, 3, 4, 5; evolution.exe produces figures 6 and 7.
