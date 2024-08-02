      
     !It finds indices and metrics of the support domain
            
      SUBROUTINE MLS_SD_ind(idir,indse,dimSD,posk,dsk,eulin_lyn)

      USE param,  ONLY: n2, n3, n2m, n3m,  &
                        rc, zz, rm, zm
      USE mpih                        

     !indse: indices of the selected eulerian point

      implicit none
     !----------------------------------------------------------    
      integer, intent(in) :: idir, dimSD
      integer, dimension(2), intent(in) :: indse
      real, dimension(2,dimSD), intent(out) :: posk,dsk
      integer, dimension(2,dimSD), intent(out) :: eulin_lyn
     !----------------------------------------------------------
      integer :: i, j, k, ind, dimx, dimy, dimz
      integer :: ii, jj, kk, m
      integer, dimension(3,3,2) :: eulin
      real, dimension(3,3,2) :: euld, eulcor
      real, dimension(:), allocatable :: x, y, z
     !----------------------------------------------------------


     !define forcing grid 
      if (idir==1) then

        dimy=n2; dimz=n3m
        allocate(y(dimy),z(dimz))
        y=rc(1:n2)
        z=zm(1:n3m)

      elseif (idir==2) then

        dimy=n2m; dimz=n3
        allocate(y(dimy),z(dimz))  
        y=rm(1:n2m)
        z=zz(1:n3)

      elseif (idir==3) then

        dimy=n2m; dimz=n3m
        allocate(y(dimy),z(dimz))
        y=rm(1:n2m)
        z=zm(1:n3m)

      endif


     !get eulerian indices matrix
      do j=1, 3
        do k=1, 3

         !move the support domain inside the flowfield by one - Y dir 
          if (indse(1)==1) then
            jj=j+1
          elseif (indse(1)==n2m) then
            jj=j-1
          else
            jj=j
          endif

         !move the support domain inside the flowfield by one - Z dir 
          if (indse(2)==1) then
            kk=k+1
          elseif (indse(2)==n3m) then
            kk=k-1
          else
            kk=k
          endif
          
          eulin(j,k,1)=indse(1)+jj-2
          eulin(j,k,2)=indse(2)+kk-2
        
        enddo
      enddo
      

      m=0
      do j=1, 3
        do k=1, 3

          m=m+1
          
          ind=eulin(j,k,1)
          posk(1,m)=y(ind)
          dsk(1,m)=y(ind+1)-y(ind)

          ind=eulin(j,k,2)
          posk(2,m)=z(ind)
          dsk(2,m)=z(ind+1)-z(ind)

          eulin_lyn(1,m)=eulin(j,k,1)
          eulin_lyn(2,m)=eulin(j,k,2)

        enddo
      enddo

      deallocate(y,z)


      ENDSUBROUTINE
  
     !====================================================================

     !It computes the local value of the transfer operator phi for 
     !the forcing stage, hence derivatives of the shape function
     !are not evaluated
     ! SDc: support domain nodes (cell center) coordinates      
     ! SDd: support domain metrics
     ! parind: marker indices in the parametric space 

      SUBROUTINE MLS_SF(pos,parind,SDc,SDd,dimSD,philoc)

      USE mls_param
      USE param,  ONLY: wfun
      USE utils_math
      USE mpih

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: dimSD
      integer, dimension(2), intent(in) :: parind
      real, dimension(2), intent(in) :: pos   
      real, dimension (2,dimSD), intent(in) :: SDc, SDd
      real, dimension(dimSD), intent(out) :: philoc
     !----------------------------------------------------------
      integer :: info, i
      real :: phieps, sumphi, det, errsum
      real, dimension(3,3) :: gp, A
      real, dimension(3,dimSD) :: B
      real, dimension(3) :: ptA, gps
     !----------------------------------------------------------
  
      phieps=1.0e-10

     !compute MLS basis functions matrix
      CALL mls_basis_func(pos,gp)
      
     !compute matrices A and B for the transfer operator 
      CALL compute_AB(wfun,pos,dimSD,SDc,SDd,A,B)

      gps=gp(1,:)

     !compute the transfer operator phi
      CALL GEsolve(A,gps,ptA,3,info)

      philoc=matmul(ptA,B)
      
     !check partition of unity for phi
      if (info==0) then
        CALL MLS_error(pos,parind,2) 
      endif

      sumphi=sum(philoc)
      errsum=abs(1.0-abs(sumphi))
      if (errsum>phieps) then
        print*, ' ... for marker at',pos
        print*, 'partition of unity is not satisfied by', errsum
        CALL stopcase
      endif


      ENDSUBROUTINE
     
     !====================================================================
     
     !It computes the local value of the transfer operator phi for
     !the load transfer
      
      SUBROUTINE MLS_SF_d(pos,parind,SDc,SDd,dimSD,phid,errsum)
      
      USE mls_param
      USE param,  ONLY: wfun
      USE utils_math

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: dimSD
      real, intent(out) :: errsum
      real, dimension(2), intent(in) :: pos
      integer, dimension(2), intent(in) :: parind
      real, dimension (2,dimSD), intent(in) :: SDc, SDd
      real, dimension(3,dimSD), intent(out) :: phid
     !----------------------------------------------------------
      integer :: i, j, k
      real :: phieps, sumphi, dummy
      integer, dimension(3) :: info
      real, dimension(3) :: C, det, gps
      real, dimension(3,3,3) :: A
      real, dimension(3,dimSD,3) :: B
      real, dimension(3,3) :: gp
      real, dimension(3,3) :: gam, Aa
     !----------------------------------------------------------

      phieps=1.0e-10

      CALL mls_basis_func(pos,gp)

      CALL compute_AB_d(wfun,pos,dimSD,SDc,SDd,A,B)


     !compute gamma
      Aa=A(:,:,1);   gps=gp(1,:)
      CALL GEsolve(Aa,gps,gam(:,1),3,info(1)) 
     

     !compute dgamma/dy
      Aa=A(:,:,2);    C=0.0
      
      C=matmul(Aa,gam(:,1))
      C=gp(2,:)-C

      Aa=A(:,:,1)

      CALL GEsolve(Aa,C,gam(:,2),3,info(2))


     !compute dgamma/dz
      Aa=A(:,:,3);   C=0.0
      
      C=matmul(Aa,gam(:,1))
      C=gp(3,:)-C

      Aa=A(:,:,1)

      CALL GEsolve(Aa,C,gam(:,3),3,info(3))


     !compute phi and its derivatives
      do i=1, dimSD
        phid(:,i)=0.0
        do j=1, 3
          phid(1,i)=phid(1,i)+gam(j,1)*B(j,i,1)
          phid(2,i)=phid(2,i)+gam(j,2)*B(j,i,1)+gam(j,1)*B(j,i,2)  
          phid(3,i)=phid(3,i)+gam(j,3)*B(j,i,1)+gam(j,1)*B(j,i,3)  
        enddo                                                    
      enddo                                                      

     
     !checks matrix inversion and partition of unity 
      if (any(info==0)) then
        CALL MLS_error(pos,parind,2) 
      endif

      sumphi=sum(phid(1,:))
      errsum=abs(1.0-abs(sumphi))
      if (errsum>phieps) then
        print*, ' ... for probe at',pos
        print*, 'partition of unity is not satisfied by',  errsum 
        print*, ''
        CALL stopcase
      endif


      ENDSUBROUTINE
     
     !====================================================================

     !It computes linear basis functions for MLS approximation 

      SUBROUTINE mls_basis_func(pos,gp)

      implicit none
     !----------------------------------------------------------
      real, dimension(2), intent(in) :: pos
      real, dimension(3,3), intent(out) :: gp
     !----------------------------------------------------------
      integer :: i
     !---------------------------------------------------------- 

      gp=0.0
      do i=1,3
         gp(i,i)=1.0
      enddo
      gp(1,2:3)=pos  

      ENDSUBROUTINE
     
     !====================================================================

      SUBROUTINE compute_AB(wfun,pos,dimSD,SDc,SDd,A,B)

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: wfun, dimSD
      real, dimension(2), intent(in) :: pos
      real, dimension (2,dimSD), intent(in) :: SDc, SDd
      real, dimension (3,3), intent(out) :: A
      real, dimension (3,dimSD), intent(out) :: B
     !----------------------------------------------------------
      integer :: i, j, k
      real, dimension (3,dimSD) :: pk
      real, dimension (3,3) :: pp
      real, dimension (2,dimSD) :: dif
      real, dimension (dimSD) :: w
     !----------------------------------------------------------

      
     !evaluate basis function value at the support domain nodes
      do i=1, dimSD
        pk(1,i)=1.0
        pk(2:3,i)=SDc(:,i)
        dif(:,i)=pos-SDc(:,i)        
      enddo

     !compute weight function array 
      CALL weight_function(dif,SDd,w,dimSD,wfun)
    
     !compute matrix B
      do i=1,3
        do j=1, dimSD
          B(i,j)=pk(i,j)*w(j)
        enddo
      enddo

     !compute matrix A
      A=0.0

      do i=1, dimSD

        do j=1, 3
          do k=1, 3
            pp(j,k)=pk(j,i)*pk(k,i)
          enddo
        enddo
        do j=1,3
          do k=1,3
            A(j,k)=A(j,k)+w(i)*pp(j,k)
          enddo
        enddo

      enddo

      
      ENDSUBROUTINE
     
     !====================================================================

      SUBROUTINE compute_AB_d(wfun,pos,dimSD,SDc,SDd,A,B)

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: wfun, dimSD
      real, dimension(2), intent(in) :: pos
      real, dimension (2,dimSD), intent(in) :: SDc, SDd
      real, dimension (3,3,3), intent(out) :: A
      real, dimension (3,dimSD,3), intent(out) :: B
     !----------------------------------------------------------
      integer :: i, j, k, l
      real, dimension (3,dimSD) :: pk
      real, dimension (3,3) :: pp
      real, dimension (2,dimSD) :: dif
      real, dimension (dimSD,3) :: w
     !----------------------------------------------------------

      
     !evaluate basis function value at the support domain nodes
      do i=1, dimSD
        pk(1,i)=1.0
        pk(2:3,i)=SDc(:,i)
        dif(:,i)=pos-SDc(:,i)
      enddo

     !compute weight function array 
      CALL weight_function_d(dif,SDd,w,dimSD,wfun)


     !compute matrix B
      do i=1,3
        do j=1, dimSD
          do k=1, 3
            B(i,j,k)=pk(i,j)*w(j,k)
          enddo
        enddo
      enddo


     !compute matrix A
      A=0.0

      do i=1, dimSD

        do j=1, 3
          do k=1, 3
            pp(j,k)=pk(j,i)*pk(k,i)
          enddo
        enddo
        do j=1, 3
          do k=1, 3
            do l=1, 3
              A(j,k,l)=A(j,k,l)+w(i,l)*pp(j,k)
            enddo
          enddo
        enddo

      enddo

 
      ENDSUBROUTINE
     
     !====================================================================

      SUBROUTINE weight_function(dif,SDd,w,dimSD,wfun)

      USE mls_param,  ONLY: radiusSD

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: dimSD, wfun
      real, dimension (2,dimSD), intent(in) :: dif
      real, dimension (2,dimSD), intent(in) :: SDd 
      real, dimension (dimSD), intent(out) :: w
     !----------------------------------------------------------
      integer :: i, k
      real :: difc, rc, rc2, rc3, rc4
      real :: wcon
      real, dimension(3) :: wc
      real, dimension (2,dimSD) :: SDds
     !----------------------------------------------------------

      wcon=0.3    !wideness coefficient for exponential weight function

      SDds=radiusSD*SDd 
      

      do i=1, dimSD        

        wc=0.0

        do k=1, 2

          difc=dif(k,i)
          rc=abs(difc)/SDds(k,i)
          rc2=rc*rc
          rc3=rc*rc*rc
          rc4=rc*rc*rc*rc

          if (wfun==1) then      !cubic spline
          
            if (rc>1.0) then
              wc(k)=0.0
            elseif (rc<=0.5) then
              wc(k)=(2.0/3.0)-4.0*rc2+4.0*rc3
            else
              wc(k)=(4.0/3.0)-4.0*rc+4.0*rc2-(4.0/3.0)*rc3
            endif

          elseif (wfun==2) then  !exponential function

            if (rc<=1.0) then
              wc(k)=exp(-(rc/wcon)**2)
            else
              wc(k)=0.0
            endif

          elseif (wfun==3) then  !quartic spline
            
            if (rc<=1.0) then
              wc(k)=1.0-6.0*rc2+8.0*rc3-3.0*rc4
            else
              wc(k)=0.0
            endif

          elseif (wfun==4) then   !Roma function  

            if (rc>1.0) then
              wc(k)=0.0
            elseif (rc<0.5) then
              wc(k)=1.0/3.0*(1.0+sqrt(-3.0*rc2+1.0))
            else
              wc(k)=1.0/6.0*(5.0-3.0*rc-sqrt(-3.0*(1.0-rc)**2+1.0))
            endif

          endif

        enddo

        w(i)=wc(1)*wc(2)

      enddo


      ENDSUBROUTINE
     
     !==================================================================== 

     !It computes the weight function matrix containing W and its derivatives 
      
      SUBROUTINE weight_function_d(dif,SDd,w,dimSD,wfun)

      USE mls_param,  ONLY: radiusSD

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: dimSD, wfun
      real, dimension (2,dimSD), intent(in) :: dif
      real, dimension (2,dimSD), intent(in) :: SDd 
      real, dimension (dimSD,3), intent(out) :: w
     !----------------------------------------------------------
      integer :: i, k
      real :: difc, rc, rc2, rc3, rc4, eps
      real :: drdc, wcon
      real, dimension(2) :: wc, dwcdc
      real, dimension (2,dimSD) :: SDds
     !----------------------------------------------------------
      
      wcon=0.3    !wideness coefficient for exponential weight function

      SDds=radiusSD*SDd 
      
      eps=1.0e-20
      
      do i=1, dimSD        

        wc=0.0
        dwcdc=0.0

        do k=1, 2

          difc=dif(k,i)

          if(abs(difc).le.eps) then
            drdc=0.0
          else
            drdc=(difc/abs(difc))/SDds(k,i)
          endif

          rc=abs(difc)/SDds(k,i)
          rc2=rc*rc
          rc3=rc*rc*rc
          rc4=rc*rc*rc*rc

          if (wfun==1) then      !cubic spline
          
            if (rc>1.0) then
              wc(k)=0.0
              dwcdc(k)=0.0
            elseif (rc<=0.5) then
              wc(k)=(2.0/3.0)-4.0*rc2+4.0*rc3
              dwcdc(k)=(-8.0*rc+12.0*rc2)*drdc
            else
              wc(k)=(4.0/3.0)-4.0*rc+4.0*rc2-(4.0/3.0)*rc3
              dwcdc(k)=(-4.0+8.0*rc-4.0*rc2)*drdc
            endif

          elseif (wfun==2) then  !exponential bell function

            if (rc<=1.0) then
              wc(k)=exp(-(rc/wcon)**2)
              dwcdc(k)=(-1.0/(wcon**2)*2.0*rc*exp(-(rc/wcon)**2))*drdc
            else
              wc(k)=0.0
              dwcdc(k)=0.0
            endif

          elseif (wfun==3) then  !quartic spline  
            
            if (rc<=1.0) then
              wc(k)=1.0-6.0*rc2+8.0*rc3-3.0*rc4
              dwcdc(k)=-12.0*rc+24.0*rc2-12.0*rc3
            else
              wc(k)=0.0
              dwcdc(k)=0.0
            endif
          
          elseif (wfun==4) then      !Roma function
          
            if (rc>1.0) then
              wc(k)=0.0
              dwcdc(k)=0.0
            elseif (rc<0.5) then
              wc(k)=1.0/3.0*(1.0+sqrt(-3.0*rc2+1.0))
              dwcdc(k)=-rc*(-3.0*rc2-1.0)**(-0.5)
            else
              wc(k)=1.0/6.0*(5.0-3.0*rc-sqrt(-3.0*(1.0-rc)**2+1.0))
              dwcdc(k)=-(1.0-rc)/(2.0*sqrt(1.0-3.0*(1.0-rc)**2))-0.5
            endif

          endif

        enddo

        w(i,1)=wc(1)*wc(2)
        w(i,2)=wc(2)*dwcdc(1)
        w(i,3)=wc(1)*dwcdc(2)

      enddo


      ENDSUBROUTINE

     !==================================================================== 

      SUBROUTINE MLS_error(pos,parind,nn)

      USE param, ONLY: n2, n3, rc, zz
      USE mls_param
      USE utils_math

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: nn
      real, dimension(2), intent(in) :: pos
      integer, dimension(2), intent(in) :: parind
     !----------------------------------------------------------
      integer ::  indy, indz
      integer, dimension(2) :: ecell
     !----------------------------------------------------------

      CALL search_arr(rc,n2,pos(1),indy)
      CALL search_arr(zz,n3,pos(2),indz)

      ecell=[indy,indz]
      
      print*,''
      if (nn==1) then
        print*, ' Lagrangian marker indices:', parind
      else
        print*, ' probes corresponding to marker indices:', parind
      endif
      print*, ' corresponding Eulerian cell indices:', ecell
      print*, ' ...here the MLS matrix has no inverse!' 
      print*,''


      ENDSUBROUTINE

     !====================================================================


