
     !perform check and i/o operations at the end of the time step

      SUBROUTINE ttime(ns)

      USE param
      USE mls_param,  ONLY: tprtag, nbcv, nbcw
      USE mpih
      USE divergence
      USE hdf5
      USE local_arrays
      USE rbm      

      implicit none
     !----------------------------------------------------------      
      integer, intent(in) :: ns
     !----------------------------------------------------------
      integer :: i, j, k, ycut, zcut
      integer :: hdf_error, itime, tmpi
      real :: dmax, tprfi, ektot, vorti
      real :: dtime1, dtime3, tmparr(2,3)
      real, dimension(2) :: vmax, vmin
      character(10) :: nsc
      character(50) :: filename
     !----------------------------------------------------------

      
      
      CALL vmaxv(vmax)
      CALL vminv(vmin)

      CALL cfl

      cflm=cflm*dt
  
      CALL divgck(dmax)
  


      if (dmax.gt.resid) then
        if (myid.eq.0) then
          write(*,*) ' Excessive large local residue for '// &
                'mass conservation'
          write(*,*) '     divergence value =', dmax
          write(*,*)          
        endif  
        
        CALL divgloc

      endif


      if (nwrit) then
        if (mod(time,tframef).lt.dt) then
          CALL mkmov_slice(ns)
          if (myid==0) then
            write(319,786) time,cmp_np1,cmv_np1,cma_np1
  786   format(10es16.5,2x,10es16.5,2x,10es16.5,2x,10es16.5,    &
                        2x,10es16.5,2x,10es16.5,2x,10es16.5,    &
                        2x,10es16.5,2x,10es16.5,2x,10es16.5)            
            flush(319)
          endif
        endif         
      endif



     !compute line-integral vorticity 
      CALL getintvort(vorti)


      ntime2=MPI_WTIME()
      

      if (myid==0) then

        dtime3=tload(2)-tload(1)
        dtime1=tit_f(2)-tit_f(1)-fortime
        ntime=ntime2-ntime1

!        write(*,778) dtime1
        write(*,782) cflm
        write(*,783) vmax(1),vmax(2),vmin(1),vmin(2)
        write(*,784) dmax

        if (ibact) then
          write(*,785) nbcv, nbcw 
        endif
        
  778   format(3x,'Fluid step solution time [s] =',f10.5)
  779   format(3x,'Forcing time [s] =',f10.5)
  780   format(3x,'Loads computation time [s] =',f10.5)
  781   format(3x,'Step solution time [s] =',f10.5)
  782   format(3x,'max cfl =',f10.5)
  783   format(3x,'max/min vel. comp. =',f10.5,f10.5,f10.5,f10.5) 
  784   format(3x,'max divg. value = ',10es16.5)        
  785   format(3x,'forcing error L2w(q2), L2w(q3) = ',10es16.5,10es16.5)
        
      endif


      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)


     !write fluid/shell restart files
      if ( mod(time,trestart)<=dt .and. ns/=0 ) then

        if (myid==0) then
          write(*,*)
          write(*,*) ' ...writing restart files'
          write(*,*)
        endif


        CALL mpi_write_continua(ns)

        CALL writeqbns

        if (myid==0) then

          filename=trim(outfol)//'/restart/master_rigid.dat'
          open(69,file=trim(filename),action='write',status='unknown')
          write(69,*) time
          write(69,*) dt,dt_o
          write(69,*) cmp_np1
          write(69,*) cmv_np1
          write(69,*) cma_np1
          close(69)

        endif

      endif


     !compute volume-averaged kinetic energy
      CALL kinen(ektot)


     !write on fluid report and rigid body report
      if (myid==0) then

        write(115,931) time, ns, dt, cflm, dmax, ektot, vorti
 931    format(e16.9,2x,i5,2x,e16.9,2x,e16.9,2x,e16.9,2x,e16.9,2x,  &
                e16.9,2x,e16.9)

        write(39,932) time, cmp_np1, cmv_np1, cma_np1
 932    format(e16.9,2x,e16.9,2x,e16.9,2x,e16.9,2x,         &
                       e16.9,2x,e16.9,2x,e16.9,2x,         &
                       e16.9,2x,e16.9,2x,e16.9,2x)

      endif


      if (myid==0) then

        if (ibact) then

          write(*,223) cmp_np1(1),cmp_np1(2),cmp_np1(3)
          write(*,224) norm2(cmv_np1(1:3)),norm2(cma_np1(1:3))
 
 223      format(3x,'py_{cm}, pz_{cm}, pa_{cm} ='10es16.5,2x,10es16.5,2x,10es16.5)
 224      format(3x,'v_{cm}, a_{cm} ='10es16.5,2x,10es16.5)

        endif   

        write(*,781) ntime
        write(*,*)

      endif

                                                                       
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)


      ENDSUBROUTINE

     !===================================================================


