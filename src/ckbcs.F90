
     !It computes the L2 norm of the error in velocity forcing
     !after the pressure correction
      
      SUBROUTINE checkbcs

      USE mpih
      USE mls_param
      USE param,  ONLY: zz, dimSD, n2m, n3m
      USE local_arrays,  ONLY: q2, q3      
      USE mpi_param,  ONLY: kstart, kend
      USE rbm

      implicit none
     !----------------------------------------------------------
      integer :: i, k, lc, nlu, nlv, n, np
      real :: ucy, ucz
      real, dimension(dimSD) :: uk
      real, dimension(2) :: pos
      real, dimension(2,dimSD) :: SDc, SDd
      integer, dimension(dimSD) :: tag
      integer, dimension(2,dimSD) :: eul
     !----------------------------------------------------------


      bcerrv=0.0
      bcerrw=0.0

      do np=1, nlm 
            
        pos=lmcc(np,:)
          
        if (zz(kstart)<pos(2).and.pos(2)<zz(kend+1)) then
          
         !direction 1 
          eul=eulins(np,1,:,:)

          ucy=0.0
          do n=1, dimSD
            uk(n)=q2(eul(1,n),eul(2,n))
            ucy=ucy+phipred(np,1,n)*uk(n)
          enddo

          bcerrv(np)=lmvel(np,1)-ucy


         !direction 2 
          eul=eulins(np,2,:,:)

          ucz=0.0
          do n=1, dimSD
            uk(n)=q3(eul(1,n),eul(2,n))
            ucz=ucz+phipred(np,2,n)*uk(n)
          enddo

          bcerrw(np)=lmvel(np,2)-ucz

        endif

      enddo
      
      if (myid==0) then
        bcerrv_r=0.0
        bcerrw_r=0.0
      endif
 
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

     
      CALL MPI_REDUCE(bcerrv,bcerrv_r,nlm, &
                      MDP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(bcerrw,bcerrw_r,nlm, &
                      MDP,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      
      if (myid==0) then

        nbcv=norm2(bcerrv_r)/float(nlm)
        nbcw=norm2(bcerrw_r)/float(nlm)

      endif

      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)


      ENDSUBROUTINE

     !============================================================

