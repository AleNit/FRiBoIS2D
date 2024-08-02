 
     !It initializes to zero main processing arrays for 
     !fluid flow solver

      SUBROUTINE initia
      
      USE param
      USE local_arrays
      USE mpi_param,    ONLY: kstart, kend
      
      implicit none
     !----------------------------------------------------------
      integer :: j, k, kc, i
     !----------------------------------------------------------
       
      do k=kstart, kend
        do j=1, n2
          pr(j,k)=0.0
          dph(j,k)=0.0
          rhs(j,k)=0.0
          ru2(j,k)=0.0
          ru3(j,k)=0.0
          ruro(j,k)=0.0
          qcap(j,k)=0.0
          hro(j,k)=0.0
        enddo
      enddo
      
      do kc=kstart-1, kend+1
        do j=1, n2
          k = kc
          
          if(k.lt.1) k=1
          if(k.gt.n3) k=n3
          
          q2(j,k)=0.0
          q3(j,k)=0.0
          dens(j,k)=1.0
        enddo
      enddo

      do j=1, n2
        rc(j)=0.0
        rm(j)=0.0
      enddo
      do k=1, n3
        zz(k)=0.0
        zm(k)=0.0
        g3rc(k)=0.0
        g3rm(k)=0.0
      enddo
      
      ENDSUBROUTINE  

     !===================================================================
     
     !It read initial solution from input file at time 0 or after a restart
     !and broadcast data to MPI process
      
      SUBROUTINE inirea
      
      USE mpih
      USE mpi_param,  ONLY: kstart, kend, countk
      USE local_arrays
      USE outflow_vars
      USE param
      
      implicit none
     !----------------------------------------------------------
      integer :: j, i 
     !----------------------------------------------------------

      
     !read primitive variables
      CALL mpi_read_continua(n2,n3,kstart,kend,2,q2)
      CALL mpi_read_continua(n2,n3,kstart,kend,3,q3)
      CALL mpi_read_continua(n2,n3,kstart,kend,4,pr)
      
     !read convective terms (values at time n-1 are needed just for
     !Adams-Bashforth scheme) 
      CALL mpi_read_continua2(n2,n3,kstart,kend,7,ru2)
      CALL mpi_read_continua2(n2,n3,kstart,kend,8,ru3)
      

     !read inflow-outflow variables
      CALL readqbns

      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr) 

      if (myid==numtasks-1) then
        do j=1, n2
          q2(j,n3)=qb2n(j)
          q3(j,n3)=qb3n(j)
        enddo
      endif


      ENDSUBROUTINE

     !===================================================================

     !It generates initial solution 
      
      SUBROUTINE inqpr
      
      USE param
      USE local_arrays,  ONLY: q2, q3, dens, pr
      USE mpi_param,     ONLY: kstart, kend
      USE mpih
      USE outflow_vars
      
      implicit none
     !---------------------------------------------------------- 
      integer :: j, k, i
      real :: xxx, yyy, eps
      real :: velin
      integer :: kc
     !----------------------------------------------------------


      velin=1.0         !uniform inflow, this must be hard-coded

     
      qb2s=0.0; qb3s=0.0
      qb2n=0.0; qb3n=0.0
      dq2x2o=0.0; dq3x2o=0.0

      eps=1.0
      do k=kstart-1, kend+1
        kc=k
        
        if (kc.lt.1) then
          kc=1
        elseif (kc.gt.n3) then
          kc=n3
        endif

        do j=1, n2
          q2(j,k)=0.0
          q3(j,k)=velin
          pr(j,k)=0.0
        enddo
      
      enddo

      do j=1, n2m
        qb3s(j)=velin
        qb3n(j)=velin
      enddo


      ENDSUBROUTINE                                  

     !===================================================================


