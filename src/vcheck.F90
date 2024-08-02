      
     !It calculates the maximum velocity over the domain

      SUBROUTINE vmaxv(vmax)
      
      USE param
      USE local_arrays,  ONLY: q2,q3,dens
      USE mpi_param,     ONLY: kstart,kend
      USE mpih
      
      implicit none
     !---------------------------------------------------------- 
      integer :: ic, jc, kc, kp
      real    :: my_vmax2, my_vmax3 
      real, dimension(3), intent(out) :: vmax
     !----------------------------------------------------------  
      
      my_vmax2=-100.0
      my_vmax3=-100.0
      
      do kc=kstart, kend
        kp=kc+1
        do jc=1, n2m
          my_vmax2=max(my_vmax2,abs(q2(jc,kc)))
          my_vmax3=max(my_vmax3,abs(q3(jc,kc)))
        enddo      
      enddo

      
      CALL MPI_ALLREDUCE(my_vmax2,vmax(1),1,MDP,MPI_MAX, &
             MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(my_vmax3,vmax(2),1,MDP,MPI_MAX, &
             MPI_COMM_WORLD,ierr)

     
      ENDSUBROUTINE

     !===================================================================

     !It computes the minimum velocity over the domain

      SUBROUTINE vminv(vmin)
      
      USE param
      USE local_arrays,  ONLY: q2,q3,dens
      USE mpi_param,     ONLY: kstart,kend
      USE mpih
      
      implicit none
     !----------------------------------------------------------
      real :: my_vmin2, my_vmin3
      integer :: jc, kc, kp, ic
      real, dimension(2), intent(out) :: vmin
     !----------------------------------------------------------

      my_vmin2=100.d0
      my_vmin3=100.d0

      do kc=kstart, kend
        kp=kc+1
        do jc=1, n2m
          my_vmin2=min(my_vmin2,q2(jc,kc))
          my_vmin3=min(my_vmin3,q3(jc,kc))
        enddo
      enddo

      CALL MPI_ALLREDUCE(my_vmin2,vmin(1),1,MDP,MPI_MIN, &
                         MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(my_vmin3,vmin(2),1,MDP,MPI_MIN, &
                         MPI_COMM_WORLD,ierr)


      ENDSUBROUTINE     

     !===================================================================


