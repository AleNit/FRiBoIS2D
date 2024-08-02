
     !It defines some metric quantities on a staggered grid 
     !   rc, zz: coordinates of the node in dir Y, Z
     !   rm, zm: coordinates of the face-center in dir Y, Z
     
     !refinement options (istr)
     ! 0 - uniform spacing
     ! 1 - wall refinement with tanh
     ! 2 - wall refinement with cos
     ! 3 - smoothed cubic B-spline stretching
     ! 100 - read grid
            
      SUBROUTINE cordin
      
      USE param                                                 
      USE mpih
      USE mls_param, ONLY: kprobe, ddmin
      
      implicit none
     !----------------------------------------------------------
      integer  :: j, k, n3mo, n2mo, nclip, i
      integer :: sizezz
      integer :: n3rem, n3pt_u, n3pt_l
      integer :: n2rem, n2pt_u, n2pt_l
      real :: tstr3, tstr2, z2dp, y2dp 
      real :: x1, x2, x3, etain, delet
      real :: he, hwu, hwl
      real :: eps, dzl, dzu, str3u_nm1, str3l_nm1, dzu_nm1, dzl_nm1
      real :: str2u_nm1, str2l_nm1
      real :: dymin,dzmin,tmp
      real, dimension(:), allocatable :: zzco, zztmp, zzup, zzdw
      real, dimension(:), allocatable :: rcco, rctmp, rcup, rcdw
     !----------------------------------------------------------

     
     !----------------------- Y coordinates definition 
      if (irest) then

        CALL readgrid(2,'input_FSI/restart/ycoord.dat')

      else

     !uniform grid
      if (istr2==0) then
        
        do j=1, n2
          x2=real(j-1)/real(n2m)
          rc(j)=rext2*x2
        enddo


     !clustering at the wall with tanh
      elseif (istr2==1) then
        
        tstr2=tanh(str2)  
        rc(1)=0.0
        do j=2, n2
          y2dp=float(2*j-n2-1)/float(n2m)
          rc(j)=(1+tanh(str2*y2dp)/tstr2)*0.5*rext2
          if (rc(j).lt.0.0.or.rc(j).gt.rext2) then
            write(*,*)'Forza la griglia: ','rc(',j,')=',rc(j)
            CALL stopcase
          endif
        enddo


     !clustering at the wall with cos
      elseif (istr2==2) then
        
        nclip=int(str2)
        n2mo=n2+nclip+nclip
        do j=1, n2mo
          etaym(j)=cos(pi*(float(j)-0.5)/float(n2mo))
        enddo
        do j=1, n2
          etay(j)=etaym(j+nclip)
        enddo
        delet=etay(1)-etay(n2)
        etain=etay(1)
        do j=1, n2
          etay(j)=etay(j)/(0.5*delet)
        enddo
        rc(1)=0.0
        do j=2, n2m
          rc(j)=rext2*(1.0-etay(j))*0.5
        enddo
        rc(n2)=rext2

      
     !smoothed cubic B-Spline stretching 
      elseif (istr2==3) then

        print*, 'istr2=3 not available'
        stop

      elseif (istr2==100) then

        CALL readgrid(2,'input_FSI/ycoord.dat')

      endif

      endif


     !corresponding cell center 
      do j=1, n2m
        rm(j)=(rc(j)+rc(j+1))*0.5
        g2rm(j)=(rc(j+1)-rc(j))*dx2
      enddo
      do j=2, n2m
        g2rc(j)=(rc(j+1)-rc(j-1))*dx2*0.5
      enddo
      g2rc(1)=(rc(2)-rc(1))*dx2
      g2rc(n2)=(rc(n2)-rc(n2m))*dx2

      rm(n2)=rc(n2)    

     !coefficients for grid spacing 
      do j=1, n2m
        udx2m(j)=dx2/g2rm(j)
        udx2c(j)=dx2/g2rc(j)
      enddo
      udx2c(n2)=dx2/g2rc(n2)



     !------------------------ Z coordinates definition
      if (irest) then
        
        CALL readgrid(3,'input_FSI/restart/zcoord.dat')
        
      else

     !uniform grid
      if (istr3==0) then
        
        do k=1,n3
          x3=real(k-1)/real(n3m)
          etaz(k)=alx3*x3
          zz(k)=etaz(k)
        enddo


     !clustering at the wall with tanh
      elseif (istr3==1) then
        
        tstr3=tanh(str3)  
        zz(1)=0.0
        do k=2, n3
          z2dp=float(2*k-n3-1)/float(n3m)
          zz(k)=(1+tanh(str3*z2dp)/tstr3)*0.5*alx3
          if (zz(k).lt.0.0.or.zz(k).gt.alx3) then
            write(*,*)'Forza la griglia: ','zc(',k,')=',zz(k)
            CALL stopcase
          endif
        enddo


     !clustering at the wall with cos
      elseif(istr3==2) then
        
        nclip = int(str3)
        n3mo = n3+nclip+nclip
        do k=1, n3mo
          etazm(k)=+cos(pi*(float(k)-0.5)/float(n3mo))
        enddo
        do k=1, n3
          etaz(k)=etazm(k+nclip)
        enddo
        delet = etaz(1)-etaz(n3)
        etain = etaz(1)
        do k=1, n3
          etaz(k)=etaz(k)/(0.5*delet)
        enddo
        zz(1) = 0.0
        do k=2, n3m
          zz(k) = alx3*(1.-etaz(k))*0.5
        enddo
        zz(n3) = alx3

      
      elseif (istr3==3) then

        print*, 'instr3=3 not available'
        stop
       
      elseif (istr3==100) then

        CALL readgrid(3,'input_FSI/zcoord.dat')

      endif

      endif


     !corresponding cell center 
      do k=1, n3m
        zm(k)=(zz(k)+zz(k+1))*0.5
        g3rm(k)=(zz(k+1)-zz(k))*dx3
      enddo
      do k=2, n3m
        g3rc(k)=(zz(k+1)-zz(k-1))*dx3*0.5
      enddo
      g3rc(1)=(zz(2)-zz(1))*dx3
      g3rc(n3)= (zz(n3)-zz(n3m))*dx3

      zm(n3)=zz(n3)


     !coefficients for grid spacing 
      do k=1, n3m
        udx3m(k)=dx3/g3rm(k)
        udx3c(k)=dx3/g3rc(k)
      end do
      udx3c(n3)=dx3/g3rc(n3)

      
     !find minimum distance in the grid and scale
      dymin=1.0e+6
      do i=2,n2
        tmp=rc(i)-rc(i-1)
        dymin=min(dymin,tmp)
      enddo

      dzmin=1.0e+6
      do i=2,n3
        tmp=zz(i)-zz(i-1)
        dzmin=min(dzmin,tmp)
      enddo

      ddmin=min(dzmin,dymin)

      ddmin=ddmin*kprobe
     
      
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

     
      ENDSUBROUTINE       

     !============================================================

      SUBROUTINE readgrid(dir,str)

      USE param, ONLY: rc,zz,n2,n3
      USE mpih

      implicit none
     !---------------------------------------------------------- 
      integer, intent(in) :: dir
      integer :: i
      character(*) :: str
     !---------------------------------------------------------- 

      if (myid==0) then     
        open(888,file=str,status='old')
        if (dir==2) then
          do i=1, n2
            read(888,*) rc(i)
          enddo
        elseif (dir==3) then
          do i=1,n3
            read(888,*) zz(i) 
          enddo
        endif  
        close(888)
      endif      

      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if (dir==2) then
        CALL MPI_BCAST(rc,n2,MDP,0,MPI_COMM_WORLD,ierr) 
      elseif (dir==3) then
        CALL MPI_BCAST(zz,n3,MDP,0,MPI_COMM_WORLD,ierr)
      endif


      ENDSUBROUTINE

