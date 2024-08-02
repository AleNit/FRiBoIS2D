
     !It assembles the RHS of the Helmholtz problem

      SUBROUTINE helsolve(nst)

      USE param
      USE local_arrays
      USE mpi_param, ONLY: kstart, kend
      USE mls_param, ONLY: fey, fez
      USE mpih
      
      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: nst
      integer :: ic, jc, kc, kstartp
     !----------------------------------------------------------

     
     !direction 2
      do kc=kstart, kend
        do jc=2, n2m
          rhs(jc,kc)=rhs2(jc,kc)+fey(jc,kc)*dt
        enddo
      enddo

      CALL solq2j       !solve for q2 in dir Y

      CALL solq2k       !solve for q3 in dir Z

     
     !set boundary conditions 
      q2(n2,:)=0.0
      q2(1,:)=0.0

      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

 
     !direction 3 
      if (kstart==1) then
        kstartp=2
      else
        kstartp=kstart
      endif

      do kc=kstartp, kend
        do jc=1, n2m
          rhs(jc,kc)=rhs3(jc,kc)+fez(jc,kc)*dt
        enddo
      enddo

      CALL solq3j      !solve for q3 in dir Y
      
      CALL solq3k       !solve for q3 in dir Z


      ENDSUBROUTINE

