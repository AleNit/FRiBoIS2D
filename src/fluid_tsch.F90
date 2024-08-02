
     !It executes one time step for the fluid solver

      SUBROUTINE fluid_tsch(ns,str)

      USE param
      USE local_arrays
      USE mpih
      USE mpi_param,  ONLY: kstart, kend
      USE mls_param,  ONLY: cbody
      USE divergence
      USE outflow_vars
      USE fsicoup
      
      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: ns
      character(1), intent(in) :: str
      integer :: n, i
      real :: myt1, myt2, myt
      real :: dmax, tr1,tr2,tr3
     !----------------------------------------------------------

      
      tit_f(1)=MPI_WTIME()
      fortime=0.0


     !compute maximum cfl value 
      CALL cfl

     
     !control on time stepping parameters
      if (.not.idtv) then           !----------------------- control of CFL
        if (ns==0) then
          cflm=cflm*dt
        elseif (ns>0) then
          dt=cflmax/cflm
          if (dt.gt.dtmax) then
            dt=dtmax
          endif
        endif
        if (dt.lt.1.0e-8) then
          if (myid==0) then
            print*, ' Too small time step size: dt =', dt
            CALL stopcase
          endif
        endif  
      else                          !----------------------- control of time-step size
        cflm=cflm*dt
        if(cflm.gt.cfllim) then
          print*, ' Too large peak CFL for stability: CFL =', cflm 
          CALL stopcase
        endif
      endif  
          
      beta=dt/ren*0.5



     !start subiterations for nonlinear terms evalutation 
      do n=1,nsst                                                 
        
        if (nsst==3) then
          al=alm(n)
          ga=gam(n)
          ro=rom(n)
        elseif (nsst==1) then
          if (ns==0) then
            ga=1.0
            ro=0.0
          else
            ga=1.0+0.5*dt/dt_o
            ro=-0.5*dt/dt_o
          endif
          al=ga+ro
        endif

      
       !Inflow/outflow BCs updated through wave equation
        CALL boucdq

   
       !compute nonlinear terms 
        CALL hdnl2

        CALL hdnl3

#ifdef DEBUG
        CALL checkpr(1)
        CALL checkqcap
#endif
     
      
       !assemble the RHS 
        CALL invtr2

        CALL invtr3
        
        CALL update_both_ghosts(n2,q2f,kstart,kend)
        CALL update_both_ghosts(n2,q3f,kstart,kend)

#ifdef DEBUG
        CALL checkvel(2)
        CALL checkvel(3)
#endif
        
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)


       !################################################## FORCING
        if (ibact) then

          myt1=MPI_WTIME()

          if (n==1) then
            CALL lgrdistribution(ns)
          endif

          CALL forcing_MLS(n)

          myt2=MPI_WTIME()
          fortime=fortime+myt2-myt1

        endif
       !################################################## 


       !solve the Helmholtz problem
        CALL helsolve(n)

        CALL update_both_ghosts(n2,q2,kstart,kend)
        CALL update_both_ghosts(n2,q3,kstart,kend)

        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)


       !compute divergence of the non-solenoidal velocity field
        CALL divg 

       
       !solve Poisson equation for pressure
        CALL phcalc_outflow
        
        CALL update_both_ghosts(n2,dph,kstart,kend)
   

     
       !update velocity field
        CALL updvp           

        CALL update_both_ghosts(n2,q2,kstart,kend)
        CALL update_both_ghosts(n2,q3,kstart,kend)


       !compute final pressure field 
        CALL prcalc       

        CALL pr_submean

        CALL update_both_ghosts(n2,pr,kstart,kend)

        CALL boucqt


      enddo

     
     !compute norm L2 of error on boundary conditions 
      if (ibact)   CALL checkbcs

      tit_f(2)=MPI_WTIME()    


     !###################################### LOADS TRANSFER
      if (ibact) then
        if (cbody) then
          CALL loads_transfer_c(ns)
        else
!          CALL loads_transfer_o(ns)
        endif        
      endif


      if (fsic==0) then

      !screen output after one time-step advancement
       if (myid==0) then
         write(*,*) &
         '============================================================='
         write(*,777) time, ns, dt
         print*,''
         
  777    format(3x,'Time = ',f10.5,3x,'step = ',i8,3x,'dt = ',es10.4)
        
       endif
    
      endif


      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)


      ENDSUBROUTINE                          



