

      SUBROUTINE fluid_post

      USE param
      USE mpih
      USE hdf5
     

      implicit none
      integer :: hdf_error
      logical :: dir_e

     
      if (myid==0) then
        print*,'';  print*,'';  print*, ''
        write(*,*) '  Time stepping end, normal exit '
        write(*,*) '  Final runtime [s] =',    ttot(2)-ttot(1)
        write(*,'(a,f10.5 )') '   Simulation final time [s] =  ', time
        print*, ''
        print*,''
      endif

     
      CALL mpi_write_continua

      CALL h5close_f(hdf_error)


     !move output files to fluid output folder
      if (myid==0) then
     
        inquire(file='log.out', exist=dir_e)
        if (dir_e) then
          CALL system('cp log.out '//trim(outfol))
          CALL system('rm -rf log.out')
        endif

      endif

      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

      CALL mem_dealloc


      ENDSUBROUTINE

