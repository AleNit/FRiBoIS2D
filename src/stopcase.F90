      
     !It stop the simulation wherever it is called

      SUBROUTINE stopcase

      USE mpih

      implicit none
      integer :: i
      logical :: ex


      if (myid==0) then

       !close opened output files 
        do i=1, 1000
          inquire(unit=i,exist=ex)
          if (ex) then
            close(i)
          endif
        enddo
  
        endif

     !interrupt simulation
      CALL MPI_ABORT(MPI_COMM_WORLD,1,ierr)

      stop


      ENDSUBROUTINE




