# FRiBoIS2D
Fluid-Rigid Body Interaction Solver 2D. This is a Fortran based 2D solver to simulate elastically mounted rigid bodies subjected to an incompressible flow. To be used for research, didactical purposes, testing method enancements. It runs efficiently on local workstations. Preprocessing and postprocessing routines are written in Matlab.

# Developers
Alessandro Nitti, Polytechnic University of Bari (https://scholar.google.it/citations?user=lv1V6-4AAAAJ&hl=it&oi=ao)  
Jietuo Wang, Polytechnic University of Bari (https://scholar.google.it/citations?user=eRV4qEMAAAAJ&hl=it&oi=ao)  
Marco Donato de Tullio, Polytechnic University of Bari (https://scholar.google.it/citations?user=4ulk_wYAAAAJ&hl=it&oi=ao)  

# Methodology

# How to run a test

# Organization of the repository
./preproc/: contanins matlab pre-processing scripts that generate the test folder and the input files needed for the execution. It is used to design the Eulerian and Lagrangian grid.  
./postproc/: contanins matlab post-processing scripts that allow to visualize the simulation output. Pressure/velocity/vorticity fields are displayed along with the body position.  
./src/: contains the makefile and all the routines to execute the program.  
./moduleload.sh: bash script used to load the necessary modules and library for execution.  
./test_OscCyl1DOF/: contains the input files for test case simulating the vortex-induced vibrations of an elastically mounted circular cylinder.  
./test_GallRact1DOF/: contains the input files for test case simulating the rotational gallopping of a rectangular wing.  
./GallRact1DOF.mp4: tape showing the output of the test and a comparison with the reference data.  
./OscCyl1DOF.mp4: tape showing the output of the test and a comparison with the reference data.  
