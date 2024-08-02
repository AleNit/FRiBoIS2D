
     !pre-processing operations for the flow solver 

      SUBROUTINE fluid_pre(dmax)

      USE param
      USE local_arrays, ONLY: q2, q3, dens, pr
      USE divergence
      USE mpih
      USE mpi_param, ONLY: kstart, kend
      USE hdf5
      USE mls_param      
      USE outflow_vars

      implicit none
     !----------------------------------------------------------
      integer :: n, l
      integer :: hdf_error
      real :: dmax
     !----------------------------------------------------------
    
      
      CALL read_fluid_in
      

     !cells count and finite difference increments 
      n2m=n2-1
      n3m=n3-1
      n2mh=n2m/2+1
      n2mp=n2mh+1


      dx2=1.0/(rext2/real(n2m))
      dx3=1.0/(alx3/real(n3m))
      dx2q=dx2**2                                                      
      dx3q=dx3**2


     !compute Reynolds and Peclet numbers
      ren=1.0/visc

      pi=2.0*asin(1.0)                          


     !open output files
      if (myid==0) then
        CALL openfi    
      endif  

    
     !assign coefficients for time marching schemes
      if(nsst>1) then   
        gam(1)=8.0/15.0
        gam(2)=5.0/12.0
        gam(3)=3.0/4.0
        rom(1)=0.0
        rom(2)=-17.0/60.0
        rom(3)=-5.0/12.0
      elseif (nsst==1) then
        gam(1)=1.5
        gam(2)=0.0
        gam(3)=0.0
        rom(1)=-0.50
        rom(2)=0.0
        rom(3)=0.0
      endif

      do n=1, nsst
        alm(n)=gam(n)+rom(n)
      enddo

     
     !allocate grid dependent variables
      CALL alloc_grid_var 

     
     !distribute work among MPI tasks  
      CALL mpi_workdistribution

      CALL mem_alloc

     
     !initialize fluid processing arrays
      CALL initia

 
     !open .h5 files (library call)
      CALL h5open_f(hdf_error) 


     !define grid metrics and indices
      CALL indic                

      CALL cordin   

    
     !flow movie parameters initialization 
      if (nwrit)     CALL inimov

     
     !preliminary calculations for pressure solver
      CALL phini 


      cflm=0.0


     !create or read initial solution and perform divergence check
      if (.not.irest) then

       !generate initial solution 
        CALL inqpr

        CALL update_both_ghosts(n2,pr,kstart,kend)
        CALL update_both_ghosts(n2,q2,kstart,kend)
        CALL update_both_ghosts(n2,q3,kstart,kend)
        
      else

       !read initial solution 
        if (myid==0) then
          write(*,*)
          write(*,*) '  ...reading restart files'
        endif  

        CALL inirea

        CALL update_both_ghosts(n2,pr,kstart,kend)
        CALL update_both_ghosts(n2,q2,kstart,kend)
        CALL update_both_ghosts(n2,q3,kstart,kend)      
        
      endif

      CALL divgck(dmax)
      CALL divgloc


      if (.not.idtv) then
        cflm=cflm*dt
      endif

      
     !find coefficients for non uniform integration in space 
      CALL coetar


      ENDSUBROUTINE

     !=================================================================== 

