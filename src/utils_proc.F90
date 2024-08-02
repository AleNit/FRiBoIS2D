      
      SUBROUTINE read_input_par_RBM(rbden,rba,    &
               cmposi,cmroti,rbgravF,rbmcon,cbody,kk,cc,    &
               kk3,rtcoup,rrr,ipert)
      
      implicit none
     !----------------------------------------------------------
      logical, intent(out) :: cbody, rtcoup
      real, intent(out) :: rbden, rba(2)
      real, intent(out) :: cmposi(2), cmroti, rbgravF(2)
      integer, intent(out) :: rbmcon(3)
      real, intent(out) :: kk(3),cc(3),kk3(3),rrr,ipert
     !----------------------------------------------------------
      integer :: i
      character(100) :: row
     !----------------------------------------------------------

      open(87,file='input_FSI/rigid_par.in', &
           status="old",action="read")

      read(87,*)
      read(87,*) cbody
      read(87,*) cmposi
      read(87,*) cmroti
      read(87,*) ipert
      read(87,*) rbgravF
      read(87,*) rbden
      read(87,*) rba
      read(87,*) rbmcon
      read(87,*) kk
      read(87,*) kk3
      read(87,*) cc
      read(87,*) rtcoup
      read(87,*) rrr
  
      close(87)


      ENDSUBROUTINE

     !============================================================

     !It reads general simulation parameters in input 
      
      SUBROUTINE read_genin

      USE param
      USE mpih
      USE mls_param
      USE fsicoup
      
      implicit none
      integer :: i
      logical :: dir_e,dir_e2
      character(2) :: str
      character(80) :: stri
     

      open(61,file='input_FSI/global_par.in', &
           status='old',action='read')

      read(61,*)
      read(61,*) dt
      read(61,*) idtv 
      read(61,*) dtmax
      read(61,*) cflmax
      read(61,*) cfllim
      read(61,*) outfol
      read(61,*) irest
      read(61,*) trestart
     
      read(61,*);  read(61,*);  read(61,*)
      read(61,*) ibact

      read(61,*);  read(61,*);  read(61,*)
      read(61,*) kprobe
      read(61,*) prlgr
      read(61,*) tprtag
      read(61,*) forst
      read(61,*) radiusSD
      read(61,*) wfun

      close(61)


     !perform some checks
      if (myid==0) then

        if (dt<1.0e-8) then
          print*, 'too small time step size' 
          stop
          CALL MPI_Abort(MPI_COMM_WORLD,1,ierr)
        endif


        !check if the simulation must be restarted and update
        !output folder name
        i=1
        dir_e=.true.
        do while (dir_e)
          write(str,'(i2.2)') i          
          inquire(file=trim(outfol)//'_'//str, exist=dir_e)
          if (dir_e) then
            irest=.true.      
            stri=trim(outfol)//'_'//str//'/restart/'   
            inquire(file=trim(stri), exist=dir_e2)
            if (dir_e2) then
              CALL system('rm -rf input_FSI/restart/')      
              CALL system('mv -f '//trim(stri)//' ./input_FSI/')
            endif
            i=i+1
          else  
            outfol=trim(outfol)//'_'//str      
            CALL system('mkdir '//trim(outfol))
            dir_e=.false.
          endif         
        enddo


      endif


      dimSD=9

      CALL MPI_BCAST(irest,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(outfol,len(outfol),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)


      ENDSUBROUTINE

     !===================================================================

     !It reads fluid simulation parameters in input and perform some checks

      SUBROUTINE read_fluid_in

      USE param
      USE mpih
      
      implicit none

      open(33,file='input_FSI/fluid_par.in',  &
            status='old',action='read')

      read(33,*)
      read(33,*) alx3
      read(33,*) rext2
      read(33,*) n3, n2
      read(33,*) istr3
      read(33,*) str3u, str3l, n3strm, n3str, n3strp
      read(33,*) int3
      read(33,*) istr2
      read(33,*) str2u, str2l, n2strm, n2str, n2strp
      read(33,*) int2
      
      read(33,*);  read(33,*); read(33,*)  
      read(33,*) nsst
      read(33,*) nwrit
      read(33,*) tframef
      read(33,*) resid
      read(33,*) cou

      read(33,*);  read(33,*); read(33,*)  
      read(33,*) visc
      read(33,*) inslws
      read(33,*) inslwn
      read(33,*) q3uw
      read(33,*) q3lw

      close(33)


     !checks on input values 
      if (nsst/=1 .and. nsst/=3) then
        if (myid==0) then
          print*, '     WARNING: Number of substeps for '// &
                  'convection term evaluation out of bounds'
          CALL stopcase
        endif
      endif

      if ((istr3/=0.and.istr3/=1.and.istr3/=2.and.    &
           istr3/=3.and.istr3/=100).or.(istr2/=0.and. &
           istr2/=1.and.istr2/=2.and.istr2/=3.and.istr2/=100)) then
        if (myid==0) then
          print*, '     WARNING: input error for '//  &
                  'bunching law'
          CALL stopcase        
        endif
      endif


      if (istr2/=3) str2=str2u
      if (istr3/=3) str3=str3u


     !check on grid stretching paramenters for spline stretching
      if (istr2==3) then
        if ((n2str+n2strm+n2strp)/=n2) then
          if (myid==0) then
            print*,''
            print*, '...inconsistent grid stretching parameters '//  &
                    'in Y direction'
            print*,''
            stop
          endif
        endif
      endif
      if (istr3==3) then
        if ((n3str+n3strm+n3strp)/=n3) then
          if (myid==0) then
            print*,''
            print*, '...inconsistent grid stretching parameters '//  &
                    'in Z direction'
            print*,''
            stop
          endif
        endif
      endif

      
      ENDSUBROUTINE

     !====================================================================

     !It allocates grid-dependent variables 
      
      SUBROUTINE alloc_grid_var

      USE mpih
      USE param
      USE mls_param
      USE outflow_vars

      implicit none


     !grid parameters
      allocate(rc(1:n2), zz(1:n3))
      allocate(rm(1:n2), zm(1:n3))
      allocate(g2rc(1:n2), g2rm(1:n2))
      allocate(g3rc(1:n3), g3rm(1:n3))
      allocate(etaz(1:n3), etazm(1:n3+500))
      allocate(etay(1:n2), etaym(1:n2+500))
 

     !quantities for derivatives
      allocate(udx3c(1:n3), udx3m(1:n3))
      allocate(udx2c(1:n2), udx2m(1:n2))


     !grid indices
      allocate(jmv(1:n2),jpv(1:n2))
      allocate(jmc(1:n2),jpc(1:n2))
      allocate(jmhv(1:n2+1))
      allocate(kmc(1:n3), kpc(1:n3), kmv(1:n3))
      allocate(kpv(1:n3), kup(1:n3), kum(1:n3))


     !metrics coefficients
      allocate(ap2j(1:n2), ac2j(1:n2), am2j(1:n2)) 
      allocate(ap3j(1:n2), ac3j(1:n2), am3j(1:n2))
      allocate(apscj(1:n2), acscj(1:n2), amscj(1:n2)) 
      allocate(ap3ck(1:n3), ac3ck(1:n3), am3ck(1:n3))
      allocate(ap3sk(1:n3), ac3sk(1:n3), am3sk(1:n3))
      allocate(ap3ssk(1:n3), ac3ssk(1:n3), am3ssk(1:n3))


     !variables for FFTW and Poisson solver
      allocate(trigx1(3*n2/2+1))
      allocate(ak2(1:n2), ap(1:n2))
      allocate(amphk(1:n3), acphk(1:n3), apphk(1:n3))
      allocate(amphj(1:n2), acphj(1:n2), apphj(1:n2))


     !outflow variables
      allocate(qb2s(n2), qb2n(n2), dqb2s(n2), dqb2n(n2))
      allocate(qb3s(n2), qb3n(n2), dqb3s(n2), dqb3n(n2))
      allocate(dq2x2o(n2), dq3x2o(n2))

     !strong coupling storage variables
      allocate(qb2s_o(n2), qb2n_o(n2))
      allocate(qb3s_o(n2), qb3n_o(n2))
      allocate(dq2x2o_o(n2), dq3x2o_o(n2))


     !other variables
      allocate(denbs(1:n2), denbn(1:n2))
      allocate(iwnt(1:n2), iwst(1:n2))
      allocate(WR(n3), zmm(n3,n3), zmmt(n3,n3))
   

      ENDSUBROUTINE

     !===================================================================

      SUBROUTINE allocate_mls_local

      USE param,  ONLY: n2, n3, dimSD, ibact
      USE mls_param
      USE mpih
      USE mpi_param
      USE rbm, ONLY: nlm

      implicit none
     !----------------------------------------------------------
      integer :: merr
     !---------------------------------------------------------- 
      
      allocate(fey(n2,kstart-lvlhalo:kend+lvlhalo-1),stat=merr)
      if(merr .ne. 0) then  
        write(6,*) 'process',myid,' failed to allocate memory for forc'
        CALL MPI_Abort(MPI_COMM_WORLD,1,ierr) 
      endif       

      allocate(fez(n2,kstart-lvlhalo:kend+lvlhalo-1),stat=merr)
      if(merr .ne. 0) then
        write(6,*) 'process',myid,' failed to allocate memory for forc'
        CALL MPI_Abort(MPI_COMM_WORLD,1,ierr)    
      endif       

      fey=0.0
      fez=0.0

      if (ibact) then
      
      allocate(eulins(nlm,2,2,dimSD))
      allocate(cl(nlm,2))
      allocate(phi(nlm,2,dimSD))
      allocate(bcerrv(nlm))
      allocate(bcerrw(nlm))
      allocate(ftag2(nlm*dimSD));  ftag2=0
      allocate(ftag3(nlm*dimSD));  ftag3=0
      allocate(phipred(nlm,2,dimSD)); phipred=0.0

      allocate(phi_L(3,3,dimSD))
      allocate(eulin_L(3,2,dimSD))
      allocate(press(nlm))
      allocate(tau(nlm,2))

      if (myid==0) then
        allocate(bcerrv_r(nlm))
        allocate(bcerrw_r(nlm))
        allocate(ftag_g2(nlm*dimSD))
        allocate(ftag_g3(nlm*dimSD))
        allocate(press_r(nlm))
        allocate(tau_r(nlm,2))
      endif

      endif

      ENDSUBROUTINE

     !===================================================================

      SUBROUTINE local_bb(cpp,cppu,cppv,bbox)

      implicit none
     !---------------------------------------------------------- 
      integer, intent(in) :: cppu,cppv
      real, dimension(cppu,cppv,4), intent(in) :: cpp
      real, dimension(8,3), intent(inout) :: bbox
     !----------------------------------------------------------
      integer :: i, j
      real :: x_n, y_n, z_n
      real :: xmax, ymax, zmax
      real :: xmin, ymin, zmin
      real :: x_nm1max, y_nm1max, z_nm1max
      real :: x_nm1min, y_nm1min, z_nm1min
     !----------------------------------------------------------
      
      x_n=cpp(1,1,1);    x_nm1max=x_n;     x_nm1min=x_n
      y_n=cpp(1,1,2);    y_nm1max=y_n;     y_nm1min=y_n
      z_n=cpp(1,1,3);    z_nm1max=z_n;     z_nm1min=z_n
      
      do i=1, cppu
        do j=1, cppv
          
          x_n=cpp(i,j,1)
          y_n=cpp(i,j,2)
          z_n=cpp(i,j,3)
          
          xmax=max(x_n,x_nm1max);   xmin=min(x_n,x_nm1min)
          ymax=max(y_n,y_nm1max);   ymin=min(y_n,y_nm1min)
          zmax=max(z_n,z_nm1max);   zmin=min(z_n,z_nm1min)
          
          x_nm1max=xmax;  y_nm1max=ymax;  z_nm1max=zmax
          x_nm1min=xmin;  y_nm1min=ymin;  z_nm1min=zmin

        enddo
      enddo

      bbox=0.0

      bbox(1:4,1)=xmax;   bbox(5:8,1)=xmin
      bbox([1,4,5,8],2)=ymax;     bbox([2,3,6,7],2)=ymin
      bbox([1,2,5,6],3)=zmin;     bbox([3,4,7,8],3)=zmax


      ENDSUBROUTINE

     !===================================================================

      SUBROUTINE print_init(dmax)

      USE param
      USE mpih,  ONLY: numtasks

      implicit none
     !----------------------------------------------------------
      integer :: nt, i
      real, intent(in) :: dmax
      integer, dimension(8) :: det
      character(10), dimension(3) :: b
     !----------------------------------------------------------
 

     !************************************** LOG HEADER
      write(*,*); write(*,*); write(*,*)

      write(*,*) '   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(*,*) '   %                                            %'
      write(*,*) '   %     2D  FSI SOLVER; RIGID BODY MOTION      %' 
      write(*,*) '   %                                            %'
      write(*,*) '   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(*,*)     
      write(*,*)
      write(*,*)
      write(*,*) '  ...solution start'
      write(*,*)

    

     !************************************** SIMULATION REPORT
    
     !display parallel computing arrangement 
      nt=0
!$OMP PARALLEL    &
!$OMP REDUCTION(+:nt)
      nt=nt+1
!$OMP  END PARALLEL

      write(141,*)
      write(141,*) 'Parallel computing arrangement'
      write(141,'(4x,a,i5)') 'MPI tasks =', numtasks
      write(141,'(4x,a,i5)') 'OMP threads per task =', nt
      write(141,*)

     !fluid setting
      write(141,*) 'Fluid settings' 
      write(141,'(4x,a,10f10.5,10f10.5,10f10.5)') &
          'Fluid Domain size {X Y Z}: ', rext2, alx3
      write(141,'(4x,a,i10,i10,i10)') &
          'Grid resolution {X Y Z}: ', n2, n3 
      write(141,*) 'node count: ', float(n2*n3)/1.0e6, 'M'
      write(141,*) 'cell count: ', float(n2m*n3m)/1.0e6, 'M'
      write(141,'(4x,a,es9.3E2)') 'Reynolds number: ', ren
  
      if (nsst>1) then
        write(141,'(4x,a)') 'Time discretization scheme for convective '// &
          'terms: Kunge-Kutta'
      else
        write(141,'(4x,a)') 'Time discretization scheme for convective '// &
          'terms: Adams-Bashfort'
      endif


      ENDSUBROUTINE

     !===================================================================

     !It computes mass-averaged total kinetic energy 
      
      SUBROUTINE kinen(ektot)

      USE param
      USE local_arrays, ONLY: q2, q3
      USE mpih
      USE mpi_param, ONLY: kstart, kend

      implicit none
     !----------------------------------------------------------
      real, intent(out) :: ektot
     !----------------------------------------------------------
      integer :: i, j, k, ip, jp, kp
      real :: myek, vol, vel, vccu, vccv, vccw
     !----------------------------------------------------------

      
      ektot=0.0
      myek=0.0


      do k=kstart, kend
        kp=k+1
        do j=1, n2m
          jp=j+1

          vol=(zz(kp)-zz(k))*(rc(jp)-rc(j))
          vccv=(q2(j,k)+q2(jp,k))*0.5
          vccw=(q3(j,k)+q3(j,kp))*0.5
          vel=vccv**2+vccw**2
          myek=myek+vel*vol

        enddo
      enddo


      CALL MPI_ALLREDUCE(myek,ektot,1,MDP,MPI_SUM,  &
                MPI_COMM_WORLD,ierr)

      ektot=0.5*ektot/(rext2*alx3)      
      

      ENDSUBROUTINE

     !===================================================================

      SUBROUTINE read_rigid

      USE rbm
      USE mpih

      implicit none
      integer :: i

      if (myid==0) then

         CALL system('< ./input_FSI/lmdata.in wc -l > nlns.txt')
         open(461,file='./nlns.txt',status='old',action='read')
         read(461,*) nlm
         close(461)
         CALL system('rm -rf nlns.txt')

         nlm=nlm-1

         open(461,file='./input_FSI/lmdata.in',   &
            status='old',action='read')

         allocate(lmcc(nlm,2),lmA(nlm),lmnor(nlm,2),vert(nlm+1,2))

         do i=1, nlm+1
           read(461,*) vert(i,1),vert(i,2)
         enddo

         close(461)

        !compute lagrangian marker position and area
         do i=1,nlm
           lmcc(i,1)=(vert(i,1)+vert(i+1,1))*0.5
           lmcc(i,2)=(vert(i,2)+vert(i+1,2))*0.5
           lmA(i)=sqrt((vert(i+1,1)-vert(i,1))**2+(vert(i+1,2)-vert(i,2))**2)
         enddo

        !create normal arrays
         CALL createnorm(nlm,vert,-1.0,lmnor)

      endif

      CALL MPI_BCAST(nlm,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if (myid/=0) then
        allocate(lmcc(nlm,2),lmA(nlm),lmnor(nlm,2),vert(nlm+1,2))
      endif
      CALL MPI_BCAST(vert,(nlm+1)*2,MDP,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lmcc,nlm*2,MDP,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lmA,nlm,MDP,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lmnor,nlm*2,MDP,0,MPI_COMM_WORLD,ierr)


     !allocate useful stuff
      allocate(lmind(nlm,2,2))
      allocate(lmvel(nlm,2))
      allocate(lmacc(nlm,2))
      allocate(Ftot(nlm,2))


      ENDSUBROUTINE

     !===================================================================

      SUBROUTINE getintvort(vir)

      USE local_arrays, ONLY: q2, q3
      USE mpih
      USE param, ONLY: n2,n2m,n3,n3m,time,zz,zm,rc,rm,outfol
      USE mpi_param, ONLY: kstart, kend

      implicit none
     !----------------------------------------------------------
      real, intent(out) :: vir
     !----------------------------------------------------------
      real :: tprfi,locz,dz,dwdy,dvdz,vi
      integer :: j,itime,indz
      integer :: myidv,tag
      character(50) :: filename
      character(5) :: inst
      real :: vortline(n2m),zatmp(n3)
     !----------------------------------------------------------

      locz=18.5

      zatmp=abs(zz-locz) 
      indz=minloc(zatmp,1)

      vi=0.0

      if (kstart<=indz.and.indz<=kend) then

       !compute out-of-plane vorticity 
        dz=(zz(indz+1)-zz(indz))
        do j=1,n2m
           dwdy=(q3(j+1,indz)-q3(j,indz))/(rc(j+1)-rc(j))
           dvdz=(q2(j,indz+1)-q2(j,indz))/dz
           vortline(j)=dwdy-dvdz
        enddo
  
       !integrate over the line 
        do j=1, n2m
          vi=vi+vortline(j)*(rc(j+1)-rc(j))
        enddo 

      endif

      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

      CALL MPI_ALLREDUCE(vi,vir,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)

      ENDSUBROUTINE

     !===================================================================


