
     !advance in time rigid body motion   

      SUBROUTINE rigid(ns,str,nk)
      
      USE utils_math
      USE param,  ONLY: dt, time, tframes, ren
      USE rbm
      USE fsicoup
      USE mls_param
      USE mpih

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: ns, nk
      character(1), intent(in) :: str
     !---------------------------------------------------------- 
      integer :: np, i, j, k, n
      integer :: config
      real, dimension(2) :: posP, velP, accP, posPo
      real, dimension(2) :: crprod1, crprod2, crprod3
      real :: rho_n(2),rho_np1(2), cmrotT(2)
      real :: accl(2), acca(2), postmp(2)
      real :: tmpcon(3,3), tmparr(2), csp(3)
      real :: forcm(2),rmat(2,2),rmati(2,2)
      character(50) :: newdir
      character(20) :: patch_ID
      real :: h_on
     !----------------------------------------------------------


      if (myid==0) then

     !external force at previous time step (weak coupling)      
      CALL int_stress(nlm,press_r,tau_r,lmA,lmnor,Ftot,forcm,momcm)


     !explicit 4th order runge-kutta
      if (rtcoup) then  !coupled system
      CALL kinem_RK4i(time,dt,forcm(1),forcm(2),momcm,rbden,rba,    &
         cmp_n(1),cmp_n(2),cmp_n(3),cmposi(1),cmposi(2),cmroti,  &
         cmv_n(1),cmv_n(2),cmv_n(3),cc(1),cc(2),cc(3),kk(1),kk(2),   &
         kk(3),kk3(1),kk3(2),kk3(3),cmp_np1(1),cmp_np1(2),           &
         cmp_np1(3),cmv_np1(1),cmv_np1(2),cmv_np1(3),cma_np1(1),     &
         cma_np1(2),cma_np1(3),rrr,rtcoup)
      else  !uncoupled dofs
       CALL kinem_RK4(time,dt,forcm(1),forcm(2),momcm,rbden,rba,    &
         cmp_n(1),cmp_n(2),cmp_n(3),cmposi(1),cmposi(2),cmroti,  &
         cmv_n(1),cmv_n(2),cmv_n(3),cc(1),cc(2),cc(3),kk(1),kk(2),   &
         kk(3),kk3(1),kk3(2),kk3(3),cmp_np1(1),cmp_np1(2),           &
         cmp_np1(3),cmv_np1(1),cmv_np1(2),cmv_np1(3),cma_np1(1),     &
         cma_np1(2),cma_np1(3))
      endif


!     !implicit midpoint rule
!      csp=[cmposi(1),cmposi(2),cmroti]
!      CALL kinem_MPR(time,dt,rbden,rba,forcm(1),forcm(2),momcm,  &
!         cc(1),cc(2),cc(3),kk(1),kk(2),kk(3),csp,           &
!         cmp_n,cmv_n,cma_n,cmp_np1,cmv_np1,cma_np1)


     !apply constraints
      tmpcon(1,:)=cmp_np1
      tmpcon(2,:)=cmv_np1
      tmpcon(3,:)=cma_np1

      do i=1,3
         if (rbmcon(i)==1) then
            tmpcon(:,i)=kinori(:,i)
         endif
      enddo
      
      cmp_np1=tmpcon(1,:)
      cmv_np1=tmpcon(2,:)
      cma_np1=tmpcon(3,:)
      

      cmtra=cmp_np1(1:2)-cmp_n(1:2)
      cmrot=cmp_np1(3)-cmp_n(3)

      
     !bring kinematics to lagrangian markers
      CALL get_rotmat(cmrot,rmat,rmati)
     
      do i=1, nlm
           
        rho_n=lmcc(i,:)-cmp_n(1:2)

        posP=cmp_np1(1:2)+matmul(rmat,rho_n) 

        crprod1(1)=-cmv_np1(3)*rho_n(2)
        crprod1(2)=cmv_np1(3)*rho_n(1)
        velP=cmv_np1(1:2)+crprod1

        crprod2(1)=-cma_np1(3)*rho_n(2)
        crprod2(2)=cma_np1(3)*rho_n(1)
        crprod3(1)=-cmv_np1(3)*crprod1(2)
        crprod3(2)=cmv_np1(3)*crprod1(1)
        accP=cma_np1(1:2)+crprod2+crprod3

        lmcc(i,:)=posP
        lmvel(i,:)=velP
        lmacc(i,:)=accP

      enddo


     !update kinematics for next time-step        
      cmp_n=cmp_np1      
      cmv_n=cmv_np1
      cma_n=cma_np1
        
    
     !apply rotation to vertices and recompute 
     !normal arrays
      do i=1, nlm+1
        rho_n=vert(i,:)-cmp_n(1:2)
        posP=cmp_np1(1:2)+matmul(rmat,rho_n)
        vert(i,:)=posP
      enddo 

      CALL createnorm(nlm,vert,-1.0,lmnor)

     
      endif !end scalar section


      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

      CALL MPI_BCAST(vert,(nlm+1)*2,MDP,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lmcc,nlm*2,MDP,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lmvel,nlm*2,MDP,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lmacc,nlm*2,MDP,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lmnor,nlm*2,MDP,0,MPI_COMM_WORLD,ierr)


      ENDSUBROUTINE

     !============================================================ 

      SUBROUTINE get_rotmat(a,R,Rinv)

      USE utils_math, ONLY: smallmat_inv
      
      implicit none
     !----------------------------------------------------------
      real, intent(in) :: a
      real, intent(out) :: R(2,2), Rinv(2,2)
      integer :: info
     !----------------------------------------------------------

      R(1,1) = COS(a)  
      R(1,2) = -SIN(a)   
      R(2,1) = SIN(a) 
      R(2,2) = COS(a)   
      
      CALL smallmat_inv(2,R,Rinv,info) 
      

      ENDSUBROUTINE

!=====================================================================================================
      
      SUBROUTINE kinem_RK4(t,dt,Fy,Fz,M,rhom,rba,posy_n,posz_n,posa_n,   &
              cspy,cspz,cspa,vely,velz,velA,cy,cz,ca,ky,kz,ka,    &
             ky3,kz3,ka3,posy_np1,posz_np1,posA_np1,vely_np1,     &
             velz_np1,velA_np1,accy_np1,accz_np1,accA_np1)

      implicit none
     !---------------------------------------------------------- 
      real, intent(in) :: t,dt
      real, intent(in) :: Fy,Fz,M,rhom,rba(2)
      real, intent(in) :: posy_n,posz_n,posa_n
      real, intent(in) :: cspy,cspz,cspa,vely,velz,velA
      real, intent(in) :: cy,cz,ca,ky,kz,ka,ky3,kz3,ka3
      real, intent(out) :: posy_np1,posz_np1,posA_np1
      real, intent(out) :: vely_np1,velz_np1,velA_np1
      real, intent(out) :: accy_np1,accz_np1,accA_np1
     !---------------------------------------------------------- 
      real :: dy,dz,da,val
      real, dimension(6) :: xo,xn,K1,K2,K3,K4
      real, dimension(6) :: tmp1,tmp2,tmp3,tmp4,FF
      real, dimension(6,6) :: Mi
     !---------------------------------------------------------- 

      dy=posy_n-cspy
      dz=posz_n-cspz
      da=posa_n-cspa

      FF=[0.0,0.0,0.0,Fy,Fz,M]
      xo=[dy,dz,da,vely,velz,velA]
      Mi=0.0
      Mi(1,1)=1.0
      Mi(2,2)=1.0
      Mi(3,3)=1.0
      Mi(4,4)=1.0/(rba(1)*rhom)
      Mi(5,5)=1.0/(rba(1)*rhom)
      Mi(6,6)=1.0/(rba(2)*rhom)


      CALL modfn(K1,xo,FF,Mi,cy,cz,ca,ky,kz,ka,ky3,kz3,ka3)

      tmp1=xo+K1*dt/2.0    
      CALL modfn(K2,tmp1,FF,Mi,cy,cz,ca,ky,kz,ka,ky3,kz3,ka3)
      
      tmp2=xo+K2*dt/2.0    
      CALL modfn(K3,tmp2,FF,Mi,cy,cz,ca,ky,kz,ka,ky3,kz3,ka3)
      
      tmp3=xo+K3*dt    
      CALL modfn(K4,tmp3,FF,Mi,cy,cz,ca,ky,kz,ka,ky3,kz3,ka3)
      
      xn=xo+dt*(K1/6.0+K2/3.0+K3/3.0+K4/6.0)

      posy_np1=cspy+xn(1)
      posz_np1=cspz+xn(2)
      posA_np1=cspa+xn(3)
      vely_np1=xn(4)
      velz_np1=xn(5)
      velA_np1=xn(6)
     
     !retrieve acceleration
      accy_np1=(vely_np1-vely)/dt
      accz_np1=(velz_np1-velz)/dt
      accA_np1=(velA_np1-velA)/dt


      ENDSUBROUTINE

!=====================================================================================================
      
     
      SUBROUTINE kinem_RK4i(t,dt,Fy,Fz,M,rhom,rba,posy_n,posz_n,posa_n,   &
              cspy,cspz,cspa,vely,velz,velA,cy,cz,ca,ky,kz,ka,    &
             ky3,kz3,ka3,posy_np1,posz_np1,posA_np1,vely_np1,     &
             velz_np1,velA_np1,accy_np1,accz_np1,accA_np1,rr,cp)

      implicit none
     !---------------------------------------------------------- 
      logical, intent(in) :: cp
      real, intent(in) :: t,dt,rr
      real, intent(in) :: Fy,Fz,M,rhom,rba(2)
      real, intent(in) :: posy_n,posz_n,posa_n
      real, intent(in) :: cspy,cspz,cspa,vely,velz,velA
      real, intent(in) :: cy,cz,ca,ky,kz,ka,ky3,kz3,ka3
      real, intent(out) :: posy_np1,posz_np1,posA_np1
      real, intent(out) :: vely_np1,velz_np1,velA_np1
      real, intent(out) :: accy_np1,accz_np1,accA_np1
     !---------------------------------------------------------- 
      real :: dy,dz,da,val
      real, dimension(6) :: xo,xn,K1,K2,K3,K4
      real, dimension(6) :: tmp1,tmp2,tmp3,tmp4,FF
      real, dimension(6,6) :: Mi
     !---------------------------------------------------------- 

      dy=posy_n-cspy
      dz=posz_n-cspz
      da=posa_n-cspa

!      FF=[0.0,0.0,0.0,Fy-val*M/rr,Fz,M]   !rotolamento sulla sx
      FF=[0.0,0.0,0.0,Fy+M/rr,Fz,M]   !rotolamento sulla dx
      xo=[dy,dz,da,vely,velz,velA]
      Mi=0.0
      Mi(1,1)=1.0
      Mi(2,2)=1.0
      Mi(3,3)=1.0
      Mi(4,4)=1.0/((rba(2)/rr**2+rba(1))*rhom)
      Mi(5,5)=1.0/(rba(1)*rhom)
      Mi(6,6)=1.0/(rba(2)*rhom)


      CALL modfn(K1,xo,FF,Mi,cy,cz,ca,ky,kz,ka,ky3,kz3,ka3)

      tmp1=xo+K1*dt/2.0    
      CALL modfn(K2,tmp1,FF,Mi,cy,cz,ca,ky,kz,ka,ky3,kz3,ka3)
      
      tmp2=xo+K2*dt/2.0    
      CALL modfn(K3,tmp2,FF,Mi,cy,cz,ca,ky,kz,ka,ky3,kz3,ka3)
      
      tmp3=xo+K3*dt    
      CALL modfn(K4,tmp3,FF,Mi,cy,cz,ca,ky,kz,ka,ky3,kz3,ka3)
      
      xn=xo+dt*(K1/6.0+K2/3.0+K3/3.0+K4/6.0)

     !rotolamento sulla dx 
      posy_np1=cspy+xn(1)
      vely_np1=xn(4)
      accy_np1=(vely_np1-vely)/dt

      posz_np1=cspz
      velz_np1=0.0
      accz_np1=0.0

      posA_np1=cspa+xn(1)/rr
      velA_np1=vely_np1/rr
      accA_np1=accA_np1/rr


      ENDSUBROUTINE

!=====================================================================================================
      
      SUBROUTINE kinem_MPR(t,dt,rhom,rba,Fy,Fz,M,cy,cz,ca,ky,kz,ka,    &
         csp,pos_n,vel_n,acc_n,pos_np1,vel_np1,acc_np1)
      
      implicit none
     !---------------------------------------------------------- 
      real, intent(in) :: dt,t
      real, intent(in) :: Fy,Fz,M,rhom,rba(2)
      real, dimension(3), intent(in) :: csp,pos_n,vel_n,acc_n
      real, intent(in) :: cy,cz,ca,ky,kz,ka
      real, dimension(3), intent(out) :: pos_np1,vel_np1,acc_np1
     !---------------------------------------------------------- 
      real :: dy,dz,da
      real, dimension(3) :: dpos,Ii,LHS,RHS,cc,kk,FF
     !---------------------------------------------------------- 

      dy=pos_n(1)-csp(1)
      dz=pos_n(2)-csp(2)
      da=pos_n(3)-csp(3)
      
      dpos=[dy,dz,da]

      FF=[Fy,Fz,M]
      cc=[cy,cz,ca]
      kk=[ky,kz,ka]
      Ii=[rba(1)*rhom,rba(1)*rhom,rba(2)*rhom]

      LHS=Ii+0.5*cc*dt+0.25*dt**2*kk
      RHS=FF-(0.5*cc*dt+0.25*dt**2*kk)*acc_n-     &
            (cc+dt*kk)*vel_n-kk*dpos
      acc_np1=RHS/LHS
      vel_np1=vel_n+0.5*dt*(acc_n+acc_np1)
      pos_np1=pos_n+dt*vel_n+0.25*dt**2*(acc_n+acc_np1)


      ENDSUBROUTINE

!=============================================================================================

      SUBROUTINE modfn(K,tmp,FF,Mi,cy,cz,ca,ky,kz,ka,ky3,kz3,ka3)

      implicit none
     !---------------------------------------------------------- 
      real, intent(in) :: tmp(6),FF(6),Mi(6,6)
      real,intent(in) :: cy,ky,cz,kz,ca,ka,ky3,kz3,ka3
      real, intent(out) :: K(6)
     !---------------------------------------------------------- 
      real :: KL(6,6), RHS(6), KNL(6,6), tmp3(6)
     !---------------------------------------------------------- 

      tmp3=tmp**3

      KL(1,:)=[0.0,0.0,0.0,1.0,0.0,0.0]
      KL(2,:)=[0.0,0.0,0.0,0.0,1.0,0.0]
      KL(3,:)=[0.0,0.0,0.0,0.0,0.0,1.0]
      KL(4,:)=[-ky,0.0,0.0,-cy,0.0,0.0]
      KL(5,:)=[0.0,-kz,0.0,0.0,-cz,0.0]
      KL(6,:)=[0.0,0.0,-ka,0.0,0.0,-ca]

      KNL=0.0
      KNL(4,1)=-ky3
      KNL(5,2)=-kz3
      KNL(6,3)=-ka3

      RHS=matmul(KL,tmp)+matmul(KNL,tmp3)+FF
      K=matmul(Mi,RHS)

      ENDSUBROUTINE

!=============================================================================================

      SUBROUTINE int_stress(nlm,press,tau,lmA,lmnor,Ftot,forcm,momcm)

      USE param, ONLY: time,outfol,dt
      USE rbm, ONLY: cmp_n,rba,rbden,rbgravF,lmcc
      USE mls_param, ONLY: prlgr,tprtag

      implicit none
     !---------------------------------------------------------- 
      integer, intent(in) :: nlm
      real, intent(in) :: press(nlm),tau(nlm,2),lmA(nlm),lmnor(nlm,2)
      real, intent(out) :: Ftot(nlm,2),forcm(2),momcm
     !---------------------------------------------------------- 
      integer :: i,itime
      real, dimension(nlm,2) :: Fp,Fv
      real :: forgcm(2),radius(2),tmparr(2),mom
      real :: fy,fz,fyp,fzp,fyv,fzv
      real :: theta, nord(2)
      real :: tprfi
      character(50) :: filename
      character(5) :: sfr
     !---------------------------------------------------------- 


      do i=1, nlm        
        Fp(i,:)=-press(i)*lmnor(i,:)*lmA(i)
        Fv(i,:)=tau(i,:)*lmA(i)
      enddo

      Ftot=Fp+Fv

     !force resultant
      fyp=0.0
      fzp=0.0
      fyv=0.0
      fzv=0.0
      do i=1, nlm
        fyp=fyp+Fp(i,1)
        fzp=fzp+Fp(i,2)
        fyv=fyv+Fv(i,1)
        fzv=fzv+Fv(i,2)
      enddo
      fy=fyp+fyv
      fz=fzp+fzv

      forgcm=(1.0-1.0/rbden)*rba(1)*rbgravF

      forcm=[fy,fz]+forgcm

     !moment resultant; moment provided by the gravity force is excluded
      momcm=0.0
      do i=1, nlm
        radius=lmcc(i,1:2)-cmp_n(1:2)
        tmparr=Ftot(i,:)
        mom=radius(1)*tmparr(2)-radius(2)*tmparr(1)
        momcm=momcm+mom
      enddo


     !write to file global force parameters
      write(126,933) time,forcm,momcm,fyp,fzp,fyv,fzv
 933  format(e16.9,2x,e16.9,2x,e16.9,2x,e16.9,2x,    &
          e16.9,2x,e16.9,2x,e16.9,2x,e16.9)

     
     !write to file lagrangian marker cartesian coordinates, 
     !normal, pressure and velocity tau component
      if (prlgr) then
          if (mod((time+dt),tprtag).lt.dt) then

            tprfi = 1.0/tprtag
            itime=nint(time*tprfi)

            write(sfr,'(i5.5)') itime

            filename=trim(outfol)//'/IB_out/LM_'//sfr//'.out'

            open(46,file=trim(filename),status='unknown',   &
               access='sequential')         
            write(46,*) nlm, time
            do i=1, nlm
              write(46,934) lmcc(i,:),lmnor(i,:),Fp(i,:),Fv(i,:)
            enddo
            close(46)
       
 934    format(e16.9,2x,e16.9,2x,e16.9,2x,e16.9,2x,   &
              e16.9,2x,e16.9,2x,e16.9,2x,e16.9)

          endif
      endif


      ENDSUBROUTINE

!=============================================================================================

      SUBROUTINE createnorm(nlm,vert,or,norm)

      USE utils_math, ONLY: cross

      implicit none
     !---------------------------------------------------------- 
      integer, intent(in) :: nlm
      real, intent(in) :: or,vert(nlm+1,2)
      real, intent(out) :: norm(nlm,2)
     !----------------------------------------------------------
      integer :: i
      real :: tri_ver(9),v1(3),v2(3),tri_nor(3),den
     !----------------------------------------------------------

      do i=1, nlm
        tri_ver(1:3)=[0.0,vert(i,1),vert(i,2)]
        tri_ver(4:6)=[0.0,vert(i+1,1),vert(i+1,2)]
        tri_ver(7:9)=[or,(vert(i,1)+vert(i+1,1))*0.5,(vert(i,2)+vert(i+1,2))*0.5]
        v1=tri_ver(4:6)-tri_ver(1:3)
        v2=tri_ver(7:9)-tri_ver(1:3)
        CALL cross(v1,v2,tri_nor)
        den=sqrt(sum(tri_nor**2))
        tri_nor=tri_nor/den
        norm(i,:)=[tri_nor(2),tri_nor(3)]
      enddo   


      ENDSUBROUTINE

!=============================================================================================

