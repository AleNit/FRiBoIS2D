     
     !It solves tridiagonal system in j direction

      SUBROUTINE solq3j

      USE param
      USE local_arrays,  ONLY: rhs
      USE mpi_param,  ONLY: kstart, kend

      implicit none
     !---------------------------------------------------------- 
      integer :: jc,kc,ic,jpiv(n2m),info,n
      real :: betadx
      real, dimension(n2):: amjl,apjl,acjl,fjl
      real :: amjT(n2m-1),acjT(n2m),apjT(n2m-1),appj(n2-3)
      real :: acjl_b
     !---------------------------------------------------------- 

     
      betadx=beta*al

      do kc=kstart, kend

        do jc=1, n2m
          acjl_b=1.0/(1.0-betadx*ac3j(jc))
          apjl(jc)=-betadx*ap3j(jc)*acjl_b
          acjl(jc)=1.0
          amjl(jc)=-betadx*am3j(jc)*acjl_b
          fjl(jc)=rhs(jc,kc)*acjl_b
        enddo

        amjT=amjl(2:n2m)
        apjT=apjl(1:(n2m-1))
        acjT=acjl(1:n2m)
    
        CALL DGTTRF(n2m,amjT,acjT,apjT,appj,jpiv,info)

       !THIS MAY GIVE ARITHMETIC EXEPTIONS IF THE FLOW FIELD IS
       !INITIALIZED TO ZERO
        CALL DGTTRS('N',n2m,1,amjT,acjT,apjT,appj,jpiv,fjl,n2m,info)

        rhs(1:n2m,kc) = fjl(1:n2m)  
        
      enddo 

      
      ENDSUBROUTINE 

     !====================================================================

     !It solves tridiagonal system in j direction 
      
      SUBROUTINE solq2j
      
      USE param
      USE local_arrays,  ONLY: rhs
      USE mpi_param
      USE mpih

      implicit none
     !---------------------------------------------------------- 
      integer :: jc,kc,info,jpiv(n2),n
      real :: betadx
      real, dimension(n2):: amjl,apjl,acjl,fjl
      real :: amjT(n2-1),acjT(n2),apjT(n2-1),appj(n2-2)
      real :: acjl_b
     !---------------------------------------------------------- 

      betadx=beta*al
      
      apjl(1)=0.0
      acjl(1)=1.0
      amjl(1)=0.0
      apjl(n2)=0.0
      acjl(n2)=1.0
      amjl(n2)=0.0

      do kc=kstart, kend

        fjl(1) = 0.0
        do jc=2, n2m
          acjl_b=1.0/(1.0-betadx*ac2j(jc))
          apjl(jc)=-betadx*ap2j(jc)*acjl_b
          acjl(jc)=1.0
          amjl(jc)=-betadx*am2j(jc)*acjl_b
          fjl(jc)=rhs(jc,kc)*acjl_b
        enddo
        fjl(n2) = 0.0
        
        amjT=amjl(2:n2)
        apjT=apjl(1:(n2-1))
        acjT=acjl(1:n2)
    
        CALL DGTTRF(n2,amjT,acjT,apjT,appj,jpiv,info)

        CALL DGTTRS('N',n2,1,amjT,acjT,apjT,appj,jpiv,fjl,n2,info)

        rhs(1:n2,kc)=fjl(1:n2)  
        
      enddo 

      
      ENDSUBROUTINE

     !====================================================================

      SUBROUTINE solq2k
      
      USE param
      USE local_arrays,  ONLY: rhs,q2
      USE mpi_param
      USE mpih
      USE outflow_vars

      implicit none
     !---------------------------------------------------------- 
      real, dimension(n3) :: amkl,apkl,ackl,fkl
      integer :: jc,kc,info,ipkv(n3m),n
      real :: betadx
      real :: amkT(n3m-1),ackT(n3m),apkT(n3m-1),appk(n3-3)
      real :: ugkks, ugkkn
      real,allocatable :: rhst(:,:)
      real :: ackl_b
     !---------------------------------------------------------- 

      allocate(rhst(1:n3,jstart:jend))
      
      betadx=beta*al

      ugkks=1.0/(g3rm(1)*g3rc(1))*dx3q*2.0*betadx
      ugkkn=1.0/(g3rm(n3m)*g3rc(n3))*dx3q*2.0*betadx

      CALL PackZ_UnpackR(rhs(:,kstart:kend),     &
                    rhst(:,jstart:jend))

      do jc=jstart, jend

        do kc=1,n3m
          ackl_b=1.0/(1.0-ac3sk(kc)*betadx)
          amkl(kc)=-am3sk(kc)*betadx*ackl_b
          ackl(kc)=1.0
          apkl(kc)=-ap3sk(kc)*betadx*ackl_b
          fkl(kc)=rhst(kc,jc)*ackl_b
        enddo

        amkT=amkl(2:n3m)
        apkT=apkl(1:(n3m-1))
        ackT=ackl(1:n3m)

        CALL DGTTRF(n3m,amkT,ackT,apkT,appk,ipkv,info)

        fkl(1)=fkl(1)+dqb2s(jc)*ugkks
        fkl(n3m)=fkl(n3m)+dqb2n(jc)*ugkkn

        CALL DGTTRS('N',n3m,1,amkT,ackT,apkT,appk,ipkv,fkl,n3m,info)

        rhst(1:n3m,jc)=fkl(1:n3m)
        
      enddo

      CALL PackR_UnpackZ(rhst(:,jstart:jend),rhs(:,kstart:kend))

      do kc=kstart, kend
        do jc=1, n2m
          q2(jc,kc)=q2(jc,kc)+rhs(jc,kc)
        enddo
      enddo

      deallocate(rhst)


      ENDSUBROUTINE

     !====================================================================
      
      SUBROUTINE solq3k
      
      USE param
      USE local_arrays,  ONLY: q3, rhs
      USE mpi_param
      USE outflow_vars
      USE mpih

      implicit none
     !---------------------------------------------------------- 
      real, dimension(n3) :: amkl,apkl,ackl, fkl
      real :: amkT(n3-1),apkT(n3-1)
      real :: appk(n3-2)
      real :: ackT(n3)
      integer :: jc,kc,info,ic,n,i
      integer :: ipkv(n3)
      integer :: kendp, kstartp
      real :: betadx, ackl_b
      real, allocatable, dimension(:,:) :: rhst
     !----------------------------------------------------------


      allocate(rhst(1:n3,jstart:jend))


      CALL PackZ_UnpackR(rhs(:,kstart:kend),rhst(:,jstart:jend))

      betadx=beta*al

      amkl(1)=0.0
      apkl(1)=0.0
      ackl(1)=1.0
      amkl(n3)=0.0
      apkl(n3)=0.0
      ackl(n3)=1.0


      do jc=jstart, jend
          
        fkl(1)=dqb3s(jc) 
        do kc=2,n3m
          ackl_b=1.0/(1.0-ac3sk(kc)*betadx)
          amkl(kc)=-am3sk(kc)*betadx*ackl_b
          ackl(kc)=1.0
          apkl(kc)=-ap3sk(kc)*betadx*ackl_b
          fkl(kc)=rhst(kc,jc)*ackl_b
        enddo
        fkl(n3)=dqb3n(jc)  

        amkT=amkl(2:n3)
        apkT=apkl(1:(n3-1))
        ackT=ackl(1:n3)

        CALL DGTTRF(n3,amkT,ackT,apkT,appk,ipkv,info)

        CALL DGTTRS('N',n3,1,amkT,ackT,apkT,appk,ipkv,fkl,n3,info)

        rhst(1:n3m,jc)= fkl(1:n3m)
        
      enddo

      CALL PackR_UnpackZ(rhst(:,jstart:jend),rhs(:,kstart:kend)) 

      
      do kc=kstart, kend
        do jc=1, n2m
          q3(jc,kc)=q3(jc,kc)+rhs(jc,kc)
        enddo
      enddo

      if (kend.eq.n3m) then
        do jc=1, n2m
          q3(jc,n3)=q3(jc,n3)+dqb3n(jc)
        enddo
      endif

      deallocate(rhst)


      ENDSUBROUTINE

     !===================================================================

     !Thomas algorithm for tridiagonal matrix inversion 
      
      SUBROUTINE tripvmy_line(ami,aci,api,rrr,n1i,n1f,m1)

      implicit none
     !----------------------------------------------------------
      integer :: n1i,n1f,m1
      real, dimension(m1) :: ami,aci,api,rrr
      real, dimension(m1) :: q,s,fei
      real    :: fn,p
      real :: eps
      integer :: ia,ii,i,l
     !----------------------------------------------------------

     !define a minimum acceptable value to avoid underflow errors
     !with streatched grids when compiling with debug flags
      eps=1.0e-60

     !vectorized for right hand side and coefficients           
      ia = n1i + 1
      ii = n1i + n1f

     !INVERSION STARTS
      q(n1i) = -api(n1i)/aci(n1i)
      s(n1i) = -ami(n1i)/aci(n1i)
      fn = rrr(n1f)
      rrr(n1i) = rrr(n1i)/aci(n1i)

     !forward elimination sweep                                         
      do i=ia, n1f
        p =1.0/(aci(i)+ami(i)*q(i-1))
        q(i)=-api(i)*p
        s(i)=-ami(i)*s(i-1)*p

        s(i)=max(eps,s(i)) !otherwise s(i) will reach underflow

        rrr(i)=(rrr(i)-ami(i)*rrr(i-1))*p
      enddo

     !backward pass                                                     
      s(n1f)=1.0
      fei(n1f)=0.0

      do l=ia, n1f
        i=ii-l
        s(i)=s(i)+q(i)*s(i+1)
        fei(i)=rrr(i)+q(i)*fei(i+1)
      enddo

      rrr(n1f)=(fn-api(i)*fei(n1i)-ami(i)*fei(n1f-1))/  &
               (api(i)*s(n1i)+ami(i)*s(n1f-1)+aci(i))

     !backward elimination pass                                         
      do l=ia, n1f
        i=ii-l
        rrr(i)=rrr(n1f)*s(i)+fei(i)
      enddo


      ENDSUBROUTINE


