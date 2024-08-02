     
     !routines for debug operations

      SUBROUTINE checkvel(ii)

      USE mpih
      USE mpi_param,  ONLY: kstart, kend
      USE local_arrays
      USE param

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: ii
      integer :: i, j, k
      real :: mck1, mck2, mck3
      real :: cksum1, cksum2, cksum3
     !----------------------------------------------------------


      mck2=0.0;  mck3=0.0
      
      do k=kstart, kend
        do j=1, n2m
          mck2=mck2+q2(j,k)
          mck3=mck3+q3(j,k)
        enddo
      enddo

      CALL MPI_REDUCE(mck3,cksum3,1,MDP,MPI_SUM,0,  &
                MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(mck2,cksum2,1,MDP,MPI_SUM,0,  &
                MPI_COMM_WORLD,ierr)
        
      if (myid.eq.0) then
        print*,''
        print*,''
        write(*,'(a,i2)') '----------------------------- after INVTR',ii
        write(*,*) 'q2cksum= ',cksum2
        write(*,*) 'q3cksum= ',cksum3
      endif
     
      ENDSUBROUTINE

     !===================================================================

      SUBROUTINE checkpr(ii)

      USE mpih
      USE mpi_param,  ONLY: kstart, kend
      USE local_arrays
      USE param

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: ii
      integer :: i, j, k
      real :: mck, cksum
     !----------------------------------------------------------


      mck=0.0 

      do k=kstart, kend
        do j=1, n2m
          mck=mck+dph(j,k)
        enddo
      enddo

      CALL MPI_REDUCE(mck,cksum,1,MDP,MPI_SUM,0,  &
                MPI_COMM_WORLD,ierr)
       
      if (myid.eq.0) then
        print*,''
        if (ii==1) then
          write(*,*) '---------------------------- after HDNL2'
        elseif(ii==2) then
          write(*,*) '---------------------------- after forcing'
        else
          write(*,*) '---------------------------- after pressure '// &
                                'correction'
        endif
        write(*,*) 'dphcksum= ',cksum
      endif
     
      ENDSUBROUTINE

     !===================================================================

      SUBROUTINE checkqcap

      USE mpih
      USE mpi_param,  ONLY: kstart, kend
      USE local_arrays
      USE param

      implicit none
     !----------------------------------------------------------
      integer :: i, j, k
      integer :: kstartp
      real :: mck, cksum
     !----------------------------------------------------------


      mck=0.0 
      
      if (kstart.eq.1) then
        kstartp=2
      else
        kstartp=kstart
      endif
 
      do k=kstart, kend
        do j=1, n2m
          mck=mck+qcap(j,k)
        enddo
      enddo
        
      CALL MPI_REDUCE(mck,cksum,1,MDP,MPI_SUM,0,  &
                MPI_COMM_WORLD,ierr)
       
      if (myid.eq.0) then
        print*,''
        write(*,*) 'qcapcksum= ',cksum
      endif
     
      ENDSUBROUTINE

     !===================================================================

      SUBROUTINE checkfield

      USE mpih
      USE mpi_param,  ONLY: kstart, kend
      USE local_arrays
      USE param

      implicit none
     !----------------------------------------------------------
      integer :: i, j, k, ic, jc, kc
      integer :: kstartp
      real :: mck, cksum
      real :: mck1, mck2, mck3, cksum1, cksum2, cksum3
     !----------------------------------------------------------

      print*,''
      
      mck2=0.0
      mck3=0.0
      do k=kstart,kend
        do j=1,n2m
          mck2=mck2+q2(j,k)
          mck3=mck3+q3(j,k)
        enddo
      enddo
      
      CALL MPI_REDUCE(mck3,cksum3,1,MDP,MPI_SUM,0, &
              MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(mck2,cksum2,1,MDP,MPI_SUM,0, &
              MPI_COMM_WORLD,ierr)
      
      if (myid.eq.0) then
        write(*,*) 'q2cksum= ',cksum2
        write(*,*) 'q3cksum= ',cksum3
      endif

      ENDSUBROUTINE  


