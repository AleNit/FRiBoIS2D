      
     !It performs the calculation of dph
          
      subroutine phcalc_outflow

      USE param
      USE local_arrays, ONLY: dph
      USE mpi_param, ONLY: istart,iend,kstart,kend,dk
      USE mpih

      implicit none
     !---------------------------------------------------------- 
      integer :: i, j, k, ip, jp, kp, imh
      integer :: info, phpiv(n2m)
      real :: check_sum
      real, dimension(n2m) :: apph,amph
      real :: amphT(n2m-1), apphT(n2m-1),acphT(n2m)
      real :: appph(n2m-2),fphj(n2m)
      real :: acphT_b
      real, dimension(n2m,n3m) :: dphR, hmn
      real, dimension(n3m,n2m) :: hmnt, vmnR, vmn
      integer :: sctk
     !---------------------------------------------------------- 

      sctk=dk*n2m

      dphR=0.0        
      CALL MPI_ALLGATHER(dph(1:n2m,kstart:kend),sctk,MDP,     &
             dphR(1:n2m,:),sctk,MDP,MPI_COMM_WORLD,ierr) 

      CALL DGEMM('n','t',n2m,n3m,n3m,1.0,dphR,n2m,zmmt,n3,0.0,hmn,n2m)

      vmn=0.0 
      do k=kstart, kend

        do j=1, n2m
          acphT_b=1.0/(acphj(j)+wr(k))
          acphT(j) = 1.0                   !diagonal elements
          apph(j)= apphj(j)*acphT_b        !super-diagonal
          amph(j)= amphj(j)*acphT_b        !lower-diagonal
          fphj(j) = hmn(j,k)*acphT_b       !rhs
        enddo

        amphT=amph(2:n2m)
        apphT=apph(1:(n2m-1))

       !LU factorization of a real tridiagonal matrix
        CALL DGTTRF(n2m, amphT, acphT, apphT, appph, phpiv, info)
       
       !Solution of the factorized system 
        CALL DGTTRS('N', n2m, 1, amphT, acphT, apphT, appph, phpiv,  &
                  fphj, n2m, info)

        vmn(k,1:n2m)=fphj

      enddo

      vmnR=0.0        
      CALL MPI_ALLREDUCE(vmn,vmnR,n2m*n3m,MDP,MPI_SUM,MPI_COMM_WORLD,ierr) 

      CALL DGEMM('n','n',n3m,n2m,n3m,1.0,zmm,n3,vmnR,n3m,0.0,hmnt,n3m)


     !assign process-wise transpose 
      do k=kstart, kend
        do j=1, n2m           
          dph(j,k)=hmnt(k,j) 
        enddo
      enddo  

      ENDSUBROUTINE


