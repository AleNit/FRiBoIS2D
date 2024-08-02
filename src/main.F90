!
!     *************************************************************
!     **                                                         **
!     **                  RIGID-BODY FSI SOLVER                  **
!     **                                                         **
!     *************************************************************
!
!=====================================================================
!
!     # Alessandro Nitti, Polytechnic University of Bari,
!       alessandro.nitti@poliba.it
!     # Marco D. de Tullio, Polytechnic University of Bari,
!       marcodonato.detullio@poliba.it
!     # Jietuo Wang, Polytechnic University of Bari,
!       jietuo.wang@poliba.it
!
!=====================================================================
!       
!     > rigid body subjected to incompressible flow
!     > rigid body-dynamics is governed by a spring-mass-damper law
!     > fluid solution is obtained by a fractional step approach on a staggered
!       Cartesian grid
!     > rigid body motion is integrated by a 4th order explicit Runge-Kutta scheme
!     > solid-fluid interface is implemented by a direct-forcing immersed
!       boundary method with Moving-least-Square transfer function
!       
!=====================================================================


      PROGRAM FSI

      USE mpi_param
      USE mpih
      USE param, ONLY: time, dt, ttot, ibact, &
                       ntime1, dt_o, stopsim,irest

      implicit none
     !---------------------------------------------------------- 
      real :: dmax, i_area, erro
      integer :: ns, nz, nk
     !---------------------------------------------------------- 


     !MPI startup 
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)

      
     !pre-processing operations
      CALL read_genin

      CALL fluid_pre(dmax)
    
      if (ibact) then
        CALL read_rigid
      endif
      CALL allocate_mls_local

      if (ibact) then
        CALL rigid_pre(i_area)  
      endif
      
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
     
     
     !print simulation setting to video 
      if (myid==0) then
        CALL print_init(dmax)
      endif  


      if (.not.irest)   time=0.0


     !start time loop 
      ttot(1)=MPI_WTIME()

      stopsim=.false.  
      ns=0
      do while (.not.stopsim) 

        ntime1=MPI_WTIME()
            
        CALL fluid_tsch(ns,'a')

        if (ibact) then
          CALL rigid(ns,'a',0)      
        endif
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

        time=time+dt

        CALL ttime(ns)

        dt_o=dt
        ns=ns+1
        
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

      enddo

      ttot(2)=MPI_WTIME()

     
      CALL fluid_post
      
      CALL stopcase
      

      ENDPROGRAM

