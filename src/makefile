########################################## DEFINITIONS

EXEC = exe.e


#------------------------- compiler
FC = mpif90


#------------------------- optimization 
FFLAGS = -O3 -msse2 


#------------------------- more flags
#FFLAGS += -w 
FFLAGS += -cpp
#FFLAGS += -Wextra 
#FFLAGS += -I/usr/include


#------------------------- variables precision
FFLAGS += -fdefault-real-8
FFLAGS += -fdefault-double-8


#------------------------- include
INCL = -I/g100_work/IscrB_LeFDAV/hdf5/include


#------------------------- optional compiling 
#OP += -DSTATS
#OP += -DSTATS2
#OP += -DDEBUG
#OP += -DFORDEBUG
#OP += -DGRIDDEBUG
#OP += -DIMPMOTION
#OP += -fopenmp -DOPENMP


#------------------------- flags for debug at runtime
#FFLAGS += -g -O0 
#FFLAGS += -ffpe-trap=zero,overflow,underflow
#FFLAGS += -ffpe-trap=overflow
#FFLAGS += -ffpe-trap=invalid
#FFLAGS += -fcheck=all 
#FFLAGS += -fbacktrace -fbounds-check 


#------------------------- libreries
LIB = -L/g100_work/PROJECTS/spack/v0.16/install/0.16.2/linux-centos8-cascadelake/gcc-10.2.0/netlib-lapack-3.8.0-2h2lhpykt4anxlvqlfdqs2a6gomfgpug/lib64 -llapack -lopenblas
LIB += -L/g100_work/IscrB_LeFDAV/hdf5/lib -lz -lhdf5_fortran -lhdf5 -ldl


#------------------------- objects
OBJS_1 = declarations.o  stopcase.o  utils_math.o   \
         utils_proc.o  rigid_pre.o rigid.o             

OBJS_2 = mpi_routines.o  fileman.o  initial_sol.o  indic.o cordin.o  \
         phini.o  divergence.o  inimov.o  coetar.o  matrix_trans.o \
         fluid_pre.o  mkmov.o  HDNL.o  ckfield.o  cfl.o  \
         solve_fluid.o  INVTR.o  phcalc_outflow.o  updvp.o  helsolve.o  vcheck.F90  \
         ttime.o  fluid_post.o  PCHIC.o  forcing_MLS.o lgrdistribution.o  \
         MLSint.o  loads_transfer.o  ckbcs.o  iobcs.o  fluid_tsch.o          


default: executable 

executable:
	@echo
	@echo '--------------------------------------------------------------'
	@echo '  Making executable : ' $(EXEC)
	@echo '--------------------------------------------------------------'
	@echo
	@echo   
	make .e
	@echo
	make partclean
		


.e: main.F90 $(OBJS_1) $(OBJS_2)
	$(FC) $(FFLAGS) $(OP) $(INCL) -o $(EXEC) main.F90 $(OBJS_1) $(OBJS_2) $(LIB) 




#--------------------------------------------------- shell routines
PCHIC.o: PCHIC.f
	$(FC) $(FFLAGS) -c PCHIC.f

utils_math.o: utils_math.f90 PCHIC.o stopcase.o
	$(FC) $(FFLAGS) -c utils_math.f90

declarations.o: declarations.F90
	$(FC) $(FFLAGS) $(INCL) -c declarations.F90

utils_proc.o: utils_proc.F90 declarations.o    \
              stopcase.o
	$(FC) $(FFLAGS) $(INCL) $(OP) -c utils_proc.F90

rigid_pre.o: rigid_pre.F90 utils_math.o utils_proc.o  \
	    stopcase.o 
	$(FC) $(FFLAGS) -c rigid_pre.F90 

rigid.o: rigid.F90 utils_math.o utils_proc.o         
	$(FC) $(FFLAGS) -c rigid.F90




#--------------------------------------------------- fluid routines
stopcase.o: stopcase.F90
	$(FC) $(FFLAGS) $(INCL) -c stopcase.F90

mpi_routines.o: mpi_routines.F90 stopcase.o
	$(FC) $(FFLAGS) $(INCL) $(OP) -c mpi_routines.F90

lgrdistribution.o: lgrdistribution.F90
	$(FC) $(FFLAGS) $(OP) -c lgrdistribution.F90

fileman.o: fileman.F90
	$(FC) $(FFLAGS) $(INCL) $(OP) -c fileman.F90

initial_sol.o: initial_sol.F90 declarations.o
	$(FC) $(FFLAGS) $(OP) -c initial_sol.F90

indic.o: indic.F90 declarations.o stopcase.o
	$(FC) $(FFLAGS) $(OP) -c indic.F90

cordin.o: cordin.F90 declarations.o stopcase.o 
	$(FC) $(FFLAGS) $(OP) -c cordin.F90

matrix_trans.o: matrix_trans.F90 declarations.o
	$(FC) $(FFLAGS) $(OP) -c matrix_trans.F90

phini.o: phini.F90 declarations.o stopcase.o
	$(FC) $(FFLAGS) $(OP) -c phini.F90

divergence.o: divergence.F90 declarations.o stopcase.o
	$(FC) $(FFLAGS) $(OP) -c divergence.F90

coetar.o: coetar.F90 declarations.o stopcase.o
	$(FC) $(FFLAGS) $(OP) -c coetar.F90

inimov.o: inimov.F90 declarations.o stopcase.o
	$(FC) $(FFLAGS) $(INCL) $(OP) -c inimov.F90

fluid_pre.o: fluid_pre.F90 declarations.o utils_proc.o fileman.o initial_sol.o \
      indic.o cordin.o phini.o divergence.o inimov.o stopcase.o
	$(FC) $(FFLAGS) $(INCL) $(OP) -c fluid_pre.F90

HDNL.o: HDNL.F90 declarations.o stopcase.o
	$(FC) $(FFLAGS) $(OP) -c HDNL.F90

ckfield.o: ckfield.F90 declarations.o stopcase.o
	$(FC) $(FFLAGS) $(OP) -c ckfield.F90

solve_fluid.o: solve_fluid.F90 declarations.o stopcase.o matrix_trans.o
	$(FC) $(FFLAGS) $(OP) -c solve_fluid.F90

INVTR.o: INVTR.F90 declarations.o solve_fluid.o stopcase.o
	$(FC) $(FFLAGS) $(OP) -c INVTR.F90

helsolve.o: helsolve.F90 declarations.o solve_fluid.F90 
	$(FC) $(FFLAGS) $(OP) -c helsolve.F90

phcalc_outflow.o: phcalc_outflow.F90 declarations.o stopcase.o
	$(FC) $(FFLAGS) $(INCL) $(OP) -c phcalc_outflow.F90

updvp.o: updvp.F90 declarations.o stopcase.o
	$(FC) $(FFLAGS) $(OP) -c updvp.F90

cfl.o: cfl.F90 declarations.o stopcase.o
	$(FC) $(FFLAGS) $(OP) -c cfl.F90

fluid_tsch.o: fluid_tsch.F90 declarations.o HDNL.o ckfield.o INVTR.o \
      divergence.o phcalc_outflow.o updvp.o stopcase.o loads_transfer.o \
      utils_math.o lgrdistribution.o
	$(FC) $(FFLAGS) $(OP) -c fluid_tsch.F90

mkmov.o: mkmov.F90 declarations.o stopcase.o
	$(FC) $(FFLAGS) $(INCL) $(OP) -c mkmov.F90

vcheck.o: vcheck.F90 declarations.o
	$(FC) $(FFLAGS) -c vcheck.F90

ttime.o: ttime.F90 declarations.o mpi_routines.o cfl.o divergence.o \
      vcheck.o stopcase.o
	$(FC) $(FFLAGS) $(INCL) $(OP) -c ttime.F90

fluid_post.o: fluid_post.F90 declarations.o mpi_routines.o \
      stopcase.o
	$(FC) $(FFLAGS) $(INCL) $(OP) -c fluid_post.F90

MLSint.o: MLSint.F90 declarations.o utils_math.o stopcase.o
	$(FC) $(FFLAGS) -c MLSint.F90

forcing_MLS.o: forcing_MLS.F90 declarations.o MLSint.o utils_math.o \
               stopcase.o mpi_routines.o utils_proc.o 
	$(FC) $(FFLAGS) $(INCL) -c forcing_MLS.F90

ckbcs.o: ckbcs.F90 declarations.o MLSint.o utils_math.o \
              mpi_routines.o utils_proc.o
	$(FC) $(FFLAGS) $(INCL) -c ckbcs.F90              

loads_transfer.o: loads_transfer.F90 MLSint.o declarations.o \
              utils_math.o mpi_routines.o utils_proc.o stopcase.o                
	$(FC) $(FFLAGS) $(INCL) -c loads_transfer.F90

iobcs.o: iobcs.F90 declarations.o  
	$(FC) $(FFLAGS) $(INCL) -c iobcs.F90



clean:
	rm -rf *.o *.e *.mod tmp.*

partclean:
	rm *.mod *.o tmp.*




