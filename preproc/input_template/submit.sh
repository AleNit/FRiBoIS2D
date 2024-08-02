#!/bin/bash   
#SBATCH --time=24:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=10
#####SBATCH --partition=g100_usr_dbg   
#SBATCH --partition=g100_usr_prod   
#SBATCH --job-name=rtcylD
#SBATCH --err=error.err   
#SBATCH --out=log.out   
#SBATCH --account=IscrC_TJeRoP 
 
 
module purge   
module load profile/eng   
module load gcc/10.2.0   
module load szip/2.1.1--gcc--10.2.0   
module load openmpi/4.1.1--gcc--10.2.0-pmi-cuda-11.5.0 
module load openblas/0.3.12--gcc--10.2.0   
module load lapack/3.8.0--gcc--10.2.0   
 
export LD_LIBRARY_PATH=/g100_work/IscrB_LeFDAV/hdf5/lib/   
export LD_LIBRARY_PATH+=:/g100_work/PROJECTS/spack/v0.16/install/0.16.2/linux-centos8-cascadelake/gcc-10.2.0/netlib-lapack-3.8.0-2h2lhpykt4anxlvqlfdqs2a6gomfgpug/lib64/   
export LD_LIBRARY_PATH+=:/g100_work/PROJECTS/spack/v0.17/prod/0.17.1/install/0.17/linux-centos8-cascadelake/gcc-10.2.0/openblas-0.3.18-qhyo6ky3e3bugwei3iqhiamlw5qmoh4o/lib/   
 
export OMP_NUM_THREADS=1   
 
srun ../src/exe.e   
 
