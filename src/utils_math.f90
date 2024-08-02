     
     !routines for math operations

      MODULE utils_math
     
      integer :: n_gpu, n_gpv
      real :: t_lsolve_i, t_lsolve_f
      real, dimension(:), allocatable :: gp_u, wg_u, gp_v, wg_v
      
      CONTAINS
     
     !============================================================
     
     !It returns knot span where u lies; rounding range must be 
     !considered for the last knot

      FUNCTION findspan (u,vec,n)
      
      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: n
      real, intent(in) :: u
      real, dimension(n+1), intent(in) :: vec
     !---------------------------------------------------------- 
      integer :: i
      real :: eps, findspan
     !----------------------------------------------------------
      
      eps=1e-10
      
      if (abs(u-vec(n+1))<eps) then  !spacial case; last knot
        findspan=n
        return
      endif
       
      do i=1, size(vec)-1
        if (u<vec(i+1)) then
          findspan=i
          return
        endif
      enddo
      
      return

      ENDFUNCTION
      
     !============================================================

     !It computes the cross product between two arrays[3] 
      
      SUBROUTINE cross(a,b,c)
   
      implicit none
     !---------------------------------------------------------- 
      real, dimension(3), intent(in)  :: a, b
      real, dimension(3), intent(out) :: c
     !----------------------------------------------------------
      
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
       
      ENDSUBROUTINE
      
     !============================================================
      
     !It computes binomial coefficient of (cb_1,cb_2) 
      
      FUNCTION coeff_bin(cb_1,cb_2)
     
      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: cb_1, cb_2
     !----------------------------------------------------------
      real :: coeff_bin
      integer :: cb_1_f, cb_2_f, diff_f, i
     !----------------------------------------------------------   
      
      cb_1_f=1;  cb_2_f=1;  diff_f=1
      
      do i=1, cb_1;  cb_1_f=cb_1_f*i;  enddo
      
      do i=1, cb_2;  cb_2_f=cb_2_f*i;  enddo
      
      do i=1, (cb_1-cb_2);  diff_f=diff_f*i;  enddo
      
      coeff_bin=cb_1_f/(diff_f*cb_2_f)
       
      ENDFUNCTION
     
     !============================================================

     !Given the number of Gauss points, it compute 2 arrays of coordinates
     !and weights for Gauss-Legendre numerical quadrature
           
      SUBROUTINE gauss_legendre_pi_wi(n_gp,xabsc,weig)
 
      implicit none
     !---------------------------------------------------------- 
      integer, intent(in) :: n_gp
     !----------------------------------------------------------
      integer :: i, j, m
      real ::  p1, p2, p3, pp, z, z1, eps, pi
      real :: epsmin
      real, dimension(n_gp) :: xabsc, weig
     !----------------------------------------------------------
   
      eps=1.e-20 !precision of Newton-Raphson iterative method
      pi=acos(-1.0)
      m=int((n_gp+1)/2) !roots of Legendre polynomial are symmetric  
      
      do i = 1, m
        z=cos(pi*(i-0.25)/(n_gp+0.5))  !initial approximation 
        z1=0.0
        
        do while (abs(z-z1).gt.eps) !main loop of refinement    
          p1=1.0
          p2=0.0
          
          do j=1, n_gp
        !loop up the recurrence relation to get the 
        !Legendre polynomial evaluated at z
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/float(j)
          enddo
          
          pp=n_gp*(z*p1-p2)/(z*z-1.0) !Legendre polynomial derivative
          z1=z
          z=z1-p1/pp
        enddo
        
        xabsc(i)=-z
        xabsc(n_gp+1-i)=z                    !symmetry 
        weig(i)=2.0/((1.0-z*z)*pp*pp)
        weig(n_gp+1-i)=weig(i)               !symmetry
      enddo
        
      if (mod(n_gp,2).ne.0) then;  xabsc((n_gp+1)/2)=0.0;  endif


     !in case of degenerated element   
      epsmin=1.0e-30
      do i=1, n_gp
        weig(i)=max(epsmin,weig(i))
      enddo

      
      ENDSUBROUTINE
     
     !============================================================

     ! QL algorithm with implicit shifts, to determine the eigenvalues and
     ! eigenvectors of a real,symmetric, tridiagonal matrix.
     !    - d : on input diagonal elements array of length n 
     !          on output, it returns the eigenvalues
     !    - e : on input sub-diagonal elements of the tridiagonal matrix,
     !          with e(1) arbitrary; on output e is destroyed
     !    - z : input as identity matrix; on output k th column of z 
     !          returns the normalized eigenvector corresponding to d(k)
                
      SUBROUTINE eigenproblem(n,mat,z,d)

      implicit none
     !------------------------------------------------------------
      integer, intent(in) :: n
      real, dimension(n), intent(inout) :: d
      real, dimension(n,n), intent(in) :: mat
      real, dimension(n,n), intent(inout) :: z
     !------------------------------------------------------------
      integer :: i, k, l, m, iter
      real :: c, b, dd, f, g, p, r, s
      real, dimension(n) :: e
      real, dimension(n,n) :: A
     !------------------------------------------------------------

      A=mat
      d=0.0;  e=0.0

      CALL H_red(A,n,d,e)

      z=A
      
      do i=2, n;  e(i-1)=e(i);  enddo 
      e(n)=0.0

      do l=1, n
        iter=0
        
 111    do m=l, n-1    
          dd=abs(d(m))+abs(d(m+1))
          if (abs(e(m))+dd.eq.dd) then
            goto 222
          endif
        enddo 
        m=n

 222    if(m.ne.l)then
          if(iter.eq.30) then
            print*, 'too many iterations in eigenproblem'
            CALL stopcase
          endif

          iter=iter+1
          g=(d(l+1)-d(l))/(2.0*e(l))
          r=sqrt(g**2+1.0**2)
          g=d(m)-d(l)+e(l)/(g+sign(r,g))
          s=1.0
          c=1.0
          p=0.0
          
          do i=m-1, l, -1
            f=s*e(i)
            b=c*e(i)
            r=sqrt(f**2+g**2)
            e(i+1)=r
            
            if(r.eq.0.)then
              d(i+1)=d(i+1)-p
              e(m)=0.0
              goto 111
            endif

            s=f/r
            c=g/r
            g=d(i+1)-p
            r=(d(i)-g)*s+2.*c*b
            p=s*r
            d(i+1)=g+p
            g=c*r-b
            
            do k=1, n   !from eigenvectors
              f=z(k,i+1)
              z(k,i+1)=s*z(k,i)+c*f
              z(k,i)=c*z(k,i)-s*f
            enddo 
            
          enddo

          d(l)=d(l)-p
          e(l)=g
          e(m)=0.0

          goto 111

        endif
      
      enddo 

      ENDSUBROUTINE

     !============================================================

     ! Householder reduction of a real, symmetric, n by n matrix A.
     ! On output, A is replaced by the orthogonal matrix Q 
     ! effecting the transformation. d returns the diagonal elements
     ! of the tridiagonal matrix, and e the off-diagonal elements,
     ! with e(1)=0

     ! From Numerical Recipes in Fortran 90
      
      SUBROUTINE H_red(A,n,d,e)

      implicit none
     !------------------------------------------------------------
      integer, intent(in) :: n
      real, dimension(n), intent(inout) :: d, e
      real, dimension(n,n), intent(inout) :: A
     !------------------------------------------------------------
      integer :: i, j, k, l
      real :: f, g, h, hh, sca
     !------------------------------------------------------------

      do i=n, 2, -1
        l=i-1
        h=0.0
        sca=0.0
        if(l.gt.1)then
          do k=1, l
            sca=sca+abs(a(i,k))
          enddo 
          if(sca.eq.0.)then
            e(i)=A(i,l)
          else
            do k=1, l
              A(i,k)=a(i,k)/sca
              h=h+A(i,k)**2
            enddo
            f=A(i,l)
            g=-sign(sqrt(h),f)
            e(i)=sca*g
            h=h-f*g
            A(i,l)=f-g
            f=0.0
            do j=1, l
              A(j,i)=A(i,j)/h
              g=0.0
              do k=1, j
                g=g+A(j,k)*A(i,k)
              enddo 
              do k=j+1, l
                g=g+A(k,j)*A(i,k)
              enddo 
              e(j)=g/h
              f=f+e(j)*A(i,j)
            enddo 
            hh=f/(h+h)
            do j=1, l
              f=A(i,j)
              g=e(j)-hh*f
              e(j)=g
              do k=1, j
                A(j,k)=A(j,k)-f*e(k)-g*A(i,k)
              enddo
            enddo 
          endif
        else
          e(i)=A(i,l)
        endif
        d(i)=h
      enddo

      d(1)=0.0
      e(1)=0.0
      do i=1, n
        l=i-1
        if(d(i).ne.0.)then
          do j=1, l
            g=0.0
            do k=1, l
              g=g+A(i,k)*A(k,j)
            enddo
            do k=1, l
              A(k,j)=A(k,j)-g*A(k,i)
            enddo 
          enddo
        endif
        d(i)=A(i,i)
        A(i,i)=1.0
        do j=1, l
          A(i,j)=0.0
          A(j,i)=0.0
        enddo 
      enddo  

      ENDSUBROUTINE
      
     !============================================================ 

     !It computes an array multiplication: s=0.5*{b}^T[A]{b} 
      
      SUBROUTINE tri_prod(n,b,A,s)

      implicit none
     !------------------------------------------------------------
      integer, intent(in) :: n
      real, dimension(n), intent(in) :: b
      real, dimension(n,n), intent(in) :: A
      real, intent(out) :: s
     !------------------------------------------------------------
      integer :: i
      real, dimension(n) :: int_array
     !------------------------------------------------------------

      do i=1, n
         int_array(i)=dot_product(b,A(:,i))
      enddo
      
      s=0.5*dot_product(int_array,b)

      ENDSUBROUTINE

     !============================================================

     !It solves the linear system of equations {x}=[A]{b} 
     !with the algorithm chosen by the user; 3 routines from 
     !scalar Lapack libraries are called
      
      SUBROUTINE lin_solve(Ai,Bi,X,N,SID,dt)

      implicit none
     !------------------------------------------------------------
      real, intent(out) :: dt
      integer, intent(in) :: N, SID
      real, dimension(N,N), intent(in) :: Ai
      real, dimension(N), intent(in) :: Bi
      real, dimension(N), intent(out) :: X
     !------------------------------------------------------------
      integer :: i
      integer :: info
      real, dimension(N,N) :: A
      real, dimension(N) :: B, IPIV
      real, dimension(:), allocatable :: WORK
      character(1) :: trans
     !------------------------------------------------------------


      A=Ai;  B=Bi    ! need to preserve input arrays


      CALL cpu_time(t_lsolve_i)
      
      select case(SID)

      case(1)
     
     !################################################################
     !#
     !# LU decomposition with partial pivoting and row interchanges 
     !#
     !################################################################
      
        CALL DGESV(N,1,A,N,IPIV,B,N,info)
  
        if (info<0) then
          print*,''
          print*, 'WARNING: the i-th value of the resolving matrix '
          print*, 'has illegal value; i=', -info
          print*,''
          CALL stopcase
        elseif (info>0) then  
          print*,''
          print*, 'WARNING: the resolving matrix is singular'
          print*,''
          CALL stopcase
        endif  

        X=B

      case(2)

     !################################################################
     !#
     !# Diagonal pivoting factorization of A.
     !#
     !# A must be symmetric  
     !#
     !################################################################

        allocate(WORK(N))

        CALL DSYSV('U',N,1,A,N,IPIV,B,N,WORK,N,info)

        if (info<0) then
          print*,''
          print*, 'WARNING: the i-th value of the resolving matrix '
          print*, 'has illegal value; i=', -info
          print*,''
          CALL stopcase
        elseif (info>0) then
          print*,''
          print*, 'WARNING: the block diagonal matrix is'
          print*, 'exactly singular, so the solution could'
          print*, 'not be completeed'
          print*,''
          CALL stopcase
        endif

        X=B

      case(3)

     !################################################################
     !#
     !# Cholesky decomposition is used to factor A.
     !#
     !# A must be symmetric and positive-definite
     !#
     !################################################################


        CALL DPOSV('U',N,1,A,N,B,N,info)

        if (info<0) then
          print*,''
          print*, 'WARNING: the i-th value of the resolving matrix '
          print*, 'has illegal value; i=', -info
          print*,''
          CALL stopcase
        elseif (info>0) then
          print*,''
          print*, 'WARNING: the leading minor of order i =', info
          print*, 'of A is not positive definite, so the Cholesky'
          print*, 'factorization could not be completeed'
          print*,''
          CALL stopcase
        endif

        X=B

      end select

      CALL cpu_time(t_lsolve_f)

      dt=t_lsolve_f-t_lsolve_i


      ENDSUBROUTINE

     !============================================================
      
     !It returns the index of the real value num in the sorted 
     !array A, or the index of the first lower value if num
     !does not belong to A
      
      SUBROUTINE search_arr(A,n,num,ind)

      implicit none
     !------------------------------------------------------------
      integer, intent(in) :: n
      real, intent(in) :: num
      integer, intent(out) :: ind
      real, dimension(n), intent(in) :: A      
     !------------------------------------------------------------
      integer :: left, right, flag, mid
     !------------------------------------------------------------
      

      left=1
      flag=0
      ind=0
      right=n

      do while (left<=right)
        mid=ceiling((left+right)/2.0)
        if (A(mid)==num.or.(A(mid+1)>num.and.A(mid)<num)) then
          ind=mid
          flag=1
          exit
        elseif (A(mid)>num) then
          right=mid-1
        else
          left=mid+1
        endif
      enddo

      if (flag==0) then
        print*, ' ...linear search on real array failed'
        CALL stopcase
      endif
    
      
      ENDSUBROUTINE

     !============================================================

     !It finds the index of the first value in the unsorted integer
     !array A equal to the integer num
      
      SUBROUTINE search_arr_int(A,n,num,ind)

      implicit none
     !------------------------------------------------------------
      integer, intent(in) :: n, num
      integer, intent(out) :: ind
      integer, dimension(n), intent(in) :: A      
     !------------------------------------------------------------
      integer :: i, flag
     !------------------------------------------------------------

      
      flag=0
      
      do i=1, n
        if (num==A(i)) then
          ind=n
          flag=1
          exit
        endif
      enddo


      if (flag==0) then
        print*, ' ...linear search on integer array failed'
        CALL stopcase
      endif

      ENDSUBROUTINE
 
     !============================================================

     !It calls routines from the PCHIC library to compute a 
     !piecewise cubic Hermite interpolating function

      SUBROUTINE get_pchip(n,x,y,m,xx,f)

      implicit none
     !------------------------------------------------------------
      integer, intent(in) :: n, m
      integer, dimension(n), intent(in) :: x
      integer, dimension(m), intent(in) :: xx
      real, dimension(n), intent(in) :: y
      real, dimension(m), intent(out) :: f
     !------------------------------------------------------------
      integer :: ierr
      real, dimension(n) :: xre
      real, dimension(n) :: D
      real, dimension(m) :: xxre
     !------------------------------------------------------------

      xre=real(x)
      xxre=real(xx)

     !compute derivatives needed to determine a monotone piecewise
     !at the interpolation points
      CALL pchim(n,xre,y,D,1,ierr)

      if (ierr<0)   CALL stopcase


     !Evaluate a piecewise cubic Hermite function
      CALL pchfe(n,xre,y,D,1,.true.,m,xxre,f,ierr)

      if (ierr<0)   CALL stopcase
      

      ENDSUBROUTINE

     !============================================================
     
     !===========================================================
     ! Solutions to a system of linear equations A*x=b
     ! Method: Gauss elimination (with scaling and pivoting)
     ! Alex G. (November 2009)
     !-----------------------------------------------------------
     ! input ...
     ! a(n,n) - array of coefficients for matrix A
     ! b(n)   - array of the right hand coefficients b
     ! n      - number of equations (size of matrix A)
     ! output ...
     ! x(n)   - solutions
     ! coments ...
     ! the original arrays a(n,n) and b(n) will be destroyed 
     ! during the calculation
     !===========================================================
      
      SUBROUTINE GEsolve(ai,bi,x,n,info)
      
      implicit none 
     !---------------------------------------------------------- 
      integer, intent(in) :: n
      integer, intent(out) :: info
      real, intent(in) :: ai(n,n), bi(n)
      real, intent(out) :: x(n)
     !---------------------------------------------------------- 
      real :: s(n)
      real :: a(n,n), b(n)
      real :: c, pivot, store
      integer :: i, j, k, l
     !---------------------------------------------------------- 


     !need to preserve input arrays 
      a=ai; b=bi
      
      info=0
    
     !step 1: begin forward elimination
      do k=1, n-1
    
       !step 2: "scaling"
       !s(i) will have the largest element from row i 
        do i=k,n                       ! loop over rows
          s(i) = 0.0
          do j=k,n                    ! loop over elements of row i
            s(i) = max(s(i),abs(a(i,j)))
          enddo
        enddo

       !step 3: "pivoting 1" 
       !find a row with the largest pivoting element
        pivot = abs(a(k,k)/s(k))
        l = k
        do j=k+1,n
          if(abs(a(j,k)/s(j)) > pivot) then
            pivot = abs(a(j,k)/s(j))
            l = j
          endif
        enddo

       !Check if the system has a sigular matrix
        if(pivot == 0.0) then
          write(*,*) ' The matrix is sigular '
          info=0
          return
        endif

       !step 4: "pivoting 2" interchange rows k and l (if needed)
        if (l /= k) then
          do j=k,n
            store = a(k,j)
            a(k,j) = a(l,j)
            a(l,j) = store
          end do
          store = b(k)
          b(k) = b(l)
          b(l) = store
        end if

       !step 5: the elimination (after scaling and pivoting)
        do i=k+1,n
          c=a(i,k)/a(k,k)
          a(i,k) = 0.0
          b(i)=b(i)- c*b(k)
          do j=k+1,n
            a(i,j) = a(i,j)-c*a(k,j)
          end do
        end do
      end do

    
     !step 6: back substiturion 
      x(n) = b(n)/a(n,n)
      do i=n-1,1,-1
        c=0.0
        do j=i+1,n
          c= c + a(i,j)*x(j)
        end do 
        x(i) = (b(i)- c)/a(i,i)
      enddo

      info=1

     
      ENDSUBROUTINE

     !============================================================

     !============================================================
     ! Inverse matrix
     ! Method: Based on Doolittle LU factorization for Ax=b
     ! Alex G. December 2009
     !-----------------------------------------------------------
     ! input ...
     ! a(n,n) - array of coefficients for matrix A
     ! n      - dimension
     ! output ...
     ! c(n,n) - inverse matrix of A
     ! comments ...
     ! the original matrix a(n,n) will be destroyed 
     ! during the calculation
     !===========================================================
      
      SUBROUTINE mat_inv(ai,c,n)

      implicit none 
     !---------------------------------------------------------- 
      integer, intent(in) :: n
      real, intent(in) :: ai(n,n)
      real, intent(out) :: c(n,n)
     !---------------------------------------------------------- 
      real :: L(n,n), U(n,n), b(n), d(n), x(n), a(n,n)
      real :: coeff
      integer i, j, k
     !---------------------------------------------------------- 


     !to preserve input matrix
      a=ai


      L=0.0;   U=0.0;   b=0.0

     !step 1: forward elimination
      do k=1, n-1
        do i=k+1,n
          coeff=a(i,k)/a(k,k)
          L(i,k) = coeff
          do j=k+1,n
            a(i,j) = a(i,j)-coeff*a(k,j)
          enddo
        enddo
      enddo

     !Step 2: prepare L and U matrices 
     !L matrix is a matrix of the elimination coefficient
     ! + the diagonal elements are 1.0
      do i=1,n
        L(i,i) = 1.0
      enddo

     !U matrix is the upper triangular part of A
      do j=1, n
        do i=1, j
          U(i,j) = a(i,j)
        enddo
      enddo

     !Step 3: compute columns of the inverse matrix C
      do k=1,n
        b(k)=1.0
        d(1) = b(1)
       
       !Step 3a: Solve Ld=b using the forward substitution
        do i=2,n
          d(i)=b(i)
          do j=1,i-1
            d(i) = d(i) - L(i,j)*d(j)
          enddo
        enddo

       !Step 3b: Solve Ux=d using the back substitution
        x(n)=d(n)/U(n,n)
        do i = n-1,1,-1
          x(i) = d(i)
          do j=n,i+1,-1
            x(i)=x(i)-U(i,j)*x(j)
          enddo
          x(i) = x(i)/u(i,i)
        enddo

       !Step 3c: fill the solutions x(n) into column k of C
        do i=1,n
          c(i,k) = x(i)
        enddo
        b(k)=0.0

      enddo

      ENDSUBROUTINE

     
     !============================================================

     !It computes the condition number for inversion defined as:
     ! c=||A||*||A^-1||, where ||.|| stands for the L_2,1 norms 
     !of a square  matrix
      
      SUBROUTINE condnum(A,n,c,time)

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: n
      real, intent(in) :: A(n,n)
      real, intent(in) :: time
      real, intent(out) :: c
     !----------------------------------------------------------
      integer :: i
      real :: Ainv(n,n), col(n)
      real :: L21_A, L21_Ainv
     !----------------------------------------------------------


      CALL mat_inv(A,Ainv,n)

      do i=1, n
        col(i)=sqrt(dot_product(A(:,i),A(:,i)))
      enddo
      L21_A=sum(col)

      do i=1, n
        col(i)=sqrt(dot_product(Ainv(:,i),Ainv(:,i)))
      enddo
      L21_Ainv=sum(col)

      c=L21_A*L21_Ainv


      ENDSUBROUTINE

     !============================================================

     !It performs direct inversion of small matrices 2x2, 3x3 or 4x4
     !3x3 and 4x4 versions are based on the subroutines M33INV and 
     !M44INV by David G. Simpson

      SUBROUTINE smallmat_inv(n,A,B,info)

      implicit none
     !---------------------------------------------------------- 
      integer, intent(in) :: n
      integer, intent(out) :: info
      real, dimension(n,n), intent(in) :: A
      real, dimension(n,n), intent(out) :: B
     !----------------------------------------------------------  
      real :: detinv, eps
     !----------------------------------------------------------  

      eps=1.0e-20

      if (n<2.or.n>4) then
        print*, ' ...not a small matrix'
        CALL stopcase
      endif

      if (n==2) then
        
        detinv = (A(1,1)*A(2,2) - A(1,2)*A(2,1))        

        if (abs(detinv)<eps) then
          info=1
          print*, 'determinant too small:',detinv
          stop
        else
          info=0          
        endif

        detinv=1.0/detinv

        B(1,1) = +detinv * A(2,2)
        B(2,1) = -detinv * A(2,1)
        B(1,2) = -detinv * A(1,2)
        B(2,2) = +detinv * A(1,1)

      elseif(n==3) then

        detinv = (A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)     &
                - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)     &
                + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))
        
        if (abs(detinv)<eps) then
          info=1
          print*, 'determinant too small:',detinv
        else
          info=0          
        endif  

        detinv=1.0/detinv

        B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
        B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
        B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
        B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
        B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
        B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
        B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
        B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
        B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))

      elseif(n==4) then

        detinv = (A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+    &
                 A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*        &
                 (A(3,2)*A(4,3)-A(3,3)*A(4,2))) - A(1,2)*(A(2,1)*    &
                 (A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*       &
                 A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-        &
                 A(3,3)*A(4,1))) + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-    &
                 A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*        &
                 A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))-      &
                 A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+       &
                 A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*        &
                 (A(3,1)*A(4,2)-A(3,2)*A(4,1))))

        if (abs(detinv)<eps) then
          info=1
          print*, 'determinant too small:',detinv
        else
          info=0          
        endif 

        detinv=1.0/detinv

        B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+       &
                 A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*        &
                 (A(3,2)*A(4,3)-A(3,3)*A(4,2)))
        B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+       &
                 A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*        &
                 (A(3,3)*A(4,1)-A(3,1)*A(4,3)))
        B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+       &
                 A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*        &
                 (A(3,1)*A(4,2)-A(3,2)*A(4,1)))
        B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+       &
                 A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*        &
                 (A(3,2)*A(4,1)-A(3,1)*A(4,2)))
        B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+       &
                 A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*        &
                 (A(3,3)*A(4,2)-A(3,2)*A(4,3)))
        B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+       &
                 A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*        &
                 (A(3,1)*A(4,3)-A(3,3)*A(4,1)))
        B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+       &
                 A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*        &
                 (A(3,2)*A(4,1)-A(3,1)*A(4,2)))
        B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+       &
                 A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*        &
                 (A(3,1)*A(4,2)-A(3,2)*A(4,1)))
        B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+       &
                 A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*        &
                 (A(2,2)*A(4,3)-A(2,3)*A(4,2)))
        B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+       &
                 A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*        &
                 (A(2,3)*A(4,1)-A(2,1)*A(4,3)))
        B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+       &
                 A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*        &
                 (A(2,1)*A(4,2)-A(2,2)*A(4,1)))
        B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+       &
                 A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*        &
                 (A(2,2)*A(4,1)-A(2,1)*A(4,2)))
        B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+       &
                 A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*        &
                 (A(2,3)*A(3,2)-A(2,2)*A(3,3)))
        B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+       &
                 A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*        &
                 (A(2,1)*A(3,3)-A(2,3)*A(3,1)))
        B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+       &
                 A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*        &
                 (A(2,2)*A(3,1)-A(2,1)*A(3,2)))
        B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+       &
                 A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*        &
                 (A(2,1)*A(3,2)-A(2,2)*A(3,1)))

      endif


      ENDSUBROUTINE

     !============================================================
              
      ENDMODULE



