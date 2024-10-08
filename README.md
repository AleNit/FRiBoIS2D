# FRiBoIS2D
Fluid-Rigid Body Interaction Solver 2D. Fortran-based 2D solver to simulate elastically-mounted rigid bodies subjected to an incompressible flow. Used for research, didactic purposes, testing method enhancements. It runs efficiently on local workstations. Preprocessing and postprocessing routines are written in Matlab.

# Developers
Alessandro Nitti, Polytechnic University of Bari (https://scholar.google.it/citations?user=lv1V6-4AAAAJ&hl=it&oi=ao)  
Jietuo Wang, Polytechnic University of Bari (https://scholar.google.it/citations?user=eRV4qEMAAAAJ&hl=it&oi=ao)  
Marco Donato de Tullio, Polytechnic University of Bari (https://scholar.google.it/citations?user=4ulk_wYAAAAJ&hl=it&oi=ao)  

# Methodology
The governing equations and the setup of the test cases are described in FRiBoiS2D_notes.pdf. Details about the methodology and numerical implementation can be found in the following papers:  
https://doi.org/10.1016/j.jcp.2016.08.020  
https://doi.org/10.1016/j.jcp.2019.01.010  
https://doi.org/10.1016/j.cma.2020.112977  

# How to run a test
1. Create the desired Eulerian grid by running the script ./proproc/MAIN_<test_name>.m. Please note that the body is expected to remain in the region of uniform grid spacing.  
2. Insert the values of physical parameters in the created test folder. For instance, change the Reynolds number by modifying line 22 of the file <test_name>/input_FSI/fluid_par.in or the solid-to-fluid density ratio at line 6 of the file <test_name>/input_FSI/rigid_par.in.
3. Compile the ./src/ files via makefile  
4. Run the test in the <test_name> folder via mpirun or srun

# Organization of the repository
./preproc/: contains matlab pre-processing scripts that generate the test folder and the input files needed for the execution. It is used to design the Eulerian and Lagrangian grid.  
./postproc/: contains matlab post-processing scripts that allow to visualize the simulation output. Pressure/velocity/vorticity fields are displayed along with the body position.  
./src/: contains the makefile and all the routines to execute the program.  
./moduleload.sh: bash script used to load the necessary modules and library for execution.  
./test_OscCyl1DOF/: contains the input files for test case simulating the vortex-induced vibrations of an elastically mounted circular cylinder.  
./test_GallRact1DOF/: contains the input files for test case simulating the rotational gallopping of a rectangular wing.  
./GallRact1DOF.mp4: tape showing the output of the test and a comparison with the reference data.  
./OscCyl1DOF.mp4: tape showing the output of the test and a comparison with the reference data.  
