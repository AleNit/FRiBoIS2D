      
     !preliminary operations for the Poisson solver

      SUBROUTINE phini
      
      USE param
      USE mpih
      
      implicit none
     !---------------------------------------------------------- 
      integer :: ic, jc, j, k, kc, kcp, kcm, mw
      integer :: nt, fftw_info
      integer :: LDB, LDC, LDA, LDVL, LDVR, NRHS, LWMAX
      integer :: INFO, LWORK
      integer :: IPIV(n3)
      real :: WI(n3), bmmt(n3,n3)
      real, dimension(:,:), allocatable :: VL, VR
      real, dimension(:), allocatable :: worke
      character :: JOBVL, JOBVR
     !---------------------------------------------------------- 

      LDA=n3
      LDVL=1
      LDVR=n3
      LDB=n3
      LDC=n3
      LWMAX=n3*n3

      allocate(VL(LDVL,n3))
      allocate(VR(LDVR,n3))
      allocate(worke(LWMAX))

    
     !Initialize tridiag matrices
      CALL tridiag_matrices   

     
     !Define a identity matrix and solve a fictitious eigenproblem
     !just to get the optimal workspace MW
      bmmt=0.0
      zmm=0.0
      zmmt=0.0
      do k=1, n3m
        bmmt(k,k)=1.0
        zmm(k,k)=1.0
        zmmt(k,k)=1.0
      enddo
      

      LWORK=-1   !Query the optimal workspace
      
      CALL DGEEV('N','V',n3m,bmmt,LDA,WR,WI,VL,LDVL,   &
                 VR,LDVR,WORKE,LWORK,INFO)
      
      MW=min(LWMAX,int(WORKE(1)))


     !Z-grid coefficients for tridiagonal matrix
      do kc=1, n3m
        bmmt(kc,kc)=acphk(kc)
        kcp=kc+1
        if (kcp.le.n3m) bmmt(kc,kcp) = apphk(kc)
        kcm = kc-1
        if (kcm.ge.1) bmmt(kc,kcm) = amphk(kc) 
      enddo


      WI=0.0
      VL=0.0
      LWORK=MW


     !Perform the eigendecomposition of the coefficient matrix:
     !Solve the proper eigenproblem to find the right eigenvector matrix
      CALL DGEEV('n','v',n3m,bmmt,LDA,WR,WI,VL,LDVL,     &
                VR,LDVR,WORKE,LWORK,INFO)

      if (INFO.GT.0) then
         write(*,*)'The algorithm failed to compute eigenvalues'
         stop
      endif
      
      zmm(1:n3m,1:n3m)=VR(1:n3m,1:n3m)

     !solve linear system VR*x=I to find the inverse of the
     !eigenvectors matrix
      CALL DGESV(n3m,n3m,VR,LDA,IPIV,zmmt,LDB,INFO)

      if (INFO.gt.0) then
        print*, 'the solution of VR*zmmt=I failed'
        stop
      endif  

      deallocate(VL,VR,worke)

      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

      
      ENDSUBROUTINE

     !============================================================ 
      
     !tridiagonal matrix coefficients at each k and i 
     !Y and Z cartesian coordinates
     !acph_ : main diagonal array
     !apph_ : upper diagonal array
     !amph_ : lower diagonal array

      SUBROUTINE tridiag_matrices

      USE param

      implicit none
     !---------------------------------------------------------- 
      integer  :: jc, jp
      integer  :: kc, km, kp
      real :: ugmmm, a33icc, a33icp
      real :: a22icc, a22icp, ac2
     !---------------------------------------------------------- 

      do jc=1, n2m
        jp=jpv(jc)
        a22icc=jmc(jc)*dx2q/g2rc(jc)
        a22icp=jpc(jc)*dx2q/g2rc(jp)
        ac2=-(a22icc+a22icp)
        ugmmm=1.0/g2rm(jc)
        amphj(jc)=a22icc*ugmmm
        apphj(jc)=a22icp*ugmmm
        acphj(jc)=-(amphj(jc)+apphj(jc))
      enddo


      do kc=1, n3m
        km=kmv(kc)
        kp=kpv(kc)
        a33icc=kmc(kc)*dx3q/g3rc(kc)
        a33icp=kpc(kc)*dx3q/g3rc(kp)
        ugmmm=1.0/g3rm(kc)
        amphk(kc)=a33icc*ugmmm
        apphk(kc)=a33icp*ugmmm
        acphk(kc)=-(amphk(kc)+apphk(kc))
      enddo


      ENDSUBROUTINE


