
     !It computes forces over a closed body

      SUBROUTINE loads_transfer_c(ns)

      USE mls_param
      USE param
      USE mpi_param,  ONLY: kstart, kend
      USE local_arrays,  ONLY: q2, q3, pr
      USE mpih
      USE fsicoup, ONLY: fsic
      USE rbm
     
      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: ns
     !----------------------------------------------------------
      integer :: np, ku, kv, idir, n, i, j, k
      integer :: lc, nlt_pa, itime
      integer :: ind1, ind2, ind3, epuind
      real :: hpr, dpdnw, dpdne
      real :: lm_prb_veln_md, lm_prb_velt_md
      real :: ducdn, d2ucdn2, d3ucdn3, ducdnm
      real :: nory, norz, tany, tanz
      real :: uk, pk, tmpu, tmpv
      real :: epu_l, epurn1, epurn2, epurn3, epurn4
      integer, dimension(2) :: tmp, indse
      integer, dimension(3,2,dimSD) :: eulin 
      real, dimension(2) :: pos, point, nor, nord
      real, dimension(2) :: lmtan,lm_prb_vel,lm_prb_veln,lm_prb_velt
      real, dimension(2) :: probeoi
      real, dimension(7) :: ucy, ucz, pc
      real, dimension(2,dimSD) :: SDc, SDd  
      real, dimension(2,2) :: Ee 
     !----------------------------------------------------------

     
      tload(1)=MPI_WTIME()

        
      hpr=ddmin*kprobe
      
      ducdn=0.0
      ducdnm=0.0
      press=0.0
      dpdne=0.0
      tau=0.0
      lmtan=0.0
      if (myid==0) then
        press_r=0.0
        tau_r=0.0
      endif

      do n=1, nlm  !---------- loop over lm

        pos=lmcc(n,:)
        probeoi=pos+lmnor(n,:)*hpr

        if (zz(kstart)<pos(2).and.pos(2)<zz(kend+1)) then

          phi_L=0.0
          eulin_L=0
          do idir=1, 3 ![v,w,p]

            if (idir==3) then
              tmp=[lmind(n,2,1),lmind(n,1,2)]
            else
              tmp=lmind(n,idir,:)
            endif 
            
            CALL probeind(tmp,probeoi,indse,idir)

           !find indices and metrics of the support domain 
            CALL MLS_SD_ind(idir,indse,dimSD,SDc,SDd,eulin_L(idir,:,:))

           !compute transfer operator array phi 
            CALL MLS_SF_d(probeoi,indse,SDc,SDd,dimSD,  &
                    phi_L(idir,:,:),epu_l) 

          enddo

!------------------------------------------ new load computation model

          CALL getpbup(eulin_L,phi_L,ucy,ucz,pc)

          !obtain unit vector along tangent direction via 
          !normal unit vector at surface and velocity vector of probe
          lm_prb_vel(1)  = ucy(1)
          lm_prb_vel(2)  = ucz(1)
          lm_prb_veln_md = dot_product(lm_prb_vel(:),lmnor(n,:))
          lm_prb_veln(:) = lm_prb_veln_md * lmnor(n,:)
          lm_prb_velt(:) = lm_prb_vel(:) - lm_prb_veln(:)
          lm_prb_velt_md = sqrt(sum(lm_prb_velt(:)**2))
          if(lm_prb_velt_md.ne.0.0) then
            lmtan(:) = lm_prb_velt(:)/sqrt(sum(lm_prb_velt(:)**2))
          endif

          nory=lmnor(n,1)
          norz=lmnor(n,2)
          tany=lmtan(  1)
          tanz=lmtan(  2)

         !transfer pressure from probes to marker,eqt(20) in wang2019
          dpdnw=dot_product(lmacc(n,:),lmnor(n,:))
          dpdnw=-1.0*dpdnw
          dpdne=dot_product(pc(2:3),lmnor(n,:))
          press(n)=pc(1)-0.5*(dpdnw+dpdne)*hpr
         
         !friction stress: eqt(26) in wang2019hydrodynamic, JCP
         !ucy, ucz,pc all are arrays with 7 elements, including
         !their values [1], 2 compoments of 1st derivatives[2-y 3-z]
         !, 4 components of 2nd derivatives [4-yy 5-yz 6-zy 7-zz]
          ducdn  = (ucy(2)*tany*nory      + ucy(3)*tany*norz + &
                    ucz(2)*tanz*nory      + ucz(3)*tanz*norz)

          d3ucdn3= (pc(4) *nory*nory      + pc(5) *nory*norz + &
                    pc(6) *norz*nory      + pc(7) *norz*norz)* &
                    1.0/visc 

          d2ucdn2= (ucy(4)*tany*nory*nory + ucy(5)*tany*nory*norz + &
                    ucy(6)*tany*norz*nory + ucy(7)*tany*norz*norz + &
                    ucz(4)*tanz*nory*nory + ucz(5)*tanz*nory*norz + &
                    ucz(6)*tanz*norz*nory + ucz(7)*tanz*norz*norz)
          ducdnm  =ducdn-d2ucdn2*hpr-0.5*d3ucdn3*hpr**2
          tau(n,:)=visc*ducdnm*lmtan(:)


!------------------------------------------ new load computation model

        endif  
  
      enddo   !---------- end loop over lagrangian markers
      
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

      CALL MPI_REDUCE(press,press_r,nlm,MDP,  &
                MPI_SUM,0,MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(tau,tau_r,nlm*2,MDP,  &
                MPI_SUM,0,MPI_COMM_WORLD,ierr)
     
     
      tload(2)=MPI_WTIME()


      ENDSUBROUTINE
    
     !============================================================ 

      SUBROUTINE probeind(markind,probepos,ind,idir)

      USE param,  ONLY: n2, n3, n2m, n3m, rc, zz, rm, zm
      USE utils_math
      USE mpih

      implicit none
     !------------------------------------------------------------
      integer, intent(in) :: markind(2), idir
      real, dimension(2), intent(in) :: probepos
      integer, dimension(2), intent(out) :: ind
     !------------------------------------------------------------
      integer :: js, je, ks, ke
      integer, parameter :: off=5
      real :: yatmp(off*2+1), zatmp(off*2+1)
     !------------------------------------------------------------

      js=markind(1)-off;   je=markind(1)+off
      ks=markind(2)-off;   ke=markind(2)+off

     !manage special cases of probe outside the domani 

      if (idir==1) then
        yatmp=abs(rc(js:je)-probepos(1))
        ind(1)=minloc(yatmp,1)
        zatmp=abs(zm(ks:ke)-probepos(2))
        ind(2)=minloc(zatmp,1)
      elseif (idir==2) then
        yatmp=abs(rm(js:je)-probepos(1))
        ind(1)=minloc(yatmp,1)
        zatmp=abs(zz(ks:ke)-probepos(2))
        ind(2)=minloc(zatmp,1)
      elseif (idir==3) then
        yatmp=abs(rm(js:je)-probepos(1))
        ind(1)=minloc(yatmp,1)
        zatmp=abs(zm(ks:ke)-probepos(2))
        ind(2)=minloc(zatmp,1)
      endif

      ind(1)=ind(1)+js-1
      ind(2)=ind(2)+ks-1


      ENDSUBROUTINE

     !============================================================ 

      SUBROUTINE getpbup(eul,phi,ucy,ucz,pc)

      USE param,  ONLY: dimSD, n2m, n3m, dx2, g2rc, dx3, g3rc
      USE local_arrays,  ONLY: q2, q3, pr

      implicit none
     !------------------------------------------------------------
      integer, dimension(3,2,dimSD), intent(in) :: eul 
      real, dimension(3,3,dimSD), intent(in) :: phi 
      real, intent(out) :: ucy(7), ucz(7), pc(7)
     !------------------------------------------------------------
      integer :: n, i
      integer :: cyk, czk
      real :: pk(3), uk(3)
      integer, dimension(dimSD) :: tag
      real, dimension(2,dimSD) :: sDdr
     !------------------------------------------------------------
      ucy=0.0
      do n=1, dimSD
        cyk   =  eul(1,1,n)
        czk   =  eul(1,2,n)
        uk(1) =  q2(cyk,czk)
        uk(2) = (q2(cyk,czk) - q2(cyk-1,czk))*dx2/g2rc(cyk)
        uk(3) = (q2(cyk,czk) - q2(cyk,czk-1))*dx3/g3rc(czk)
        do i=1, 3
          ucy(i)  =ucy(i)  +phi(1,i,n)*uk(1)
        enddo
        do i=1, 2
          ucy(i+3)=ucy(i+3)+phi(1,2,n)*uk(1+i)
        enddo
        do i=1, 2
          ucy(i+5)=ucy(i+5)+phi(1,3,n)*uk(1+i)
        enddo
      enddo

      ucz=0.0
      do n=1, dimSD
        cyk   =  eul(2,1,n)
        czk   =  eul(2,2,n)
        uk(1) =  q3(cyk,czk)
        uk(2) = (q3(cyk,czk) - q3(cyk-1,czk))*dx2/g2rc(cyk)
        uk(3) = (q3(cyk,czk) - q3(cyk,czk-1))*dx3/g3rc(czk)
        do i=1, 3
          ucz(i)  =ucz(i)  +phi(2,i,n)*uk(1)
        enddo
        do i=1, 2
          ucz(i+3)=ucz(i+3)+phi(2,2,n)*uk(1+i)
        enddo
        do i=1, 2
          ucz(i+5)=ucz(i+5)+phi(2,3,n)*uk(1+i)
        enddo
      enddo

      pc=0.0
      do n=1, dimSD
        cyk   =  eul(3,1,n)
        czk   =  eul(3,2,n)
        pk(1) =  pr(cyk,czk)
        pk(2) = (pr(cyk,czk) - pr(cyk-1,czk))*dx2/g2rc(cyk)
        pk(3) = (pr(cyk,czk) - pr(cyk,czk-1))*dx3/g3rc(czk)
        do i=1, 3
          pc(i)  =pc(i)  +phi(3,i,n)*pk(1)
        enddo
        do i=1, 2
          pc(i+3)=pc(i+3)+phi(3,2,n)*pk(1+i)
        enddo
        do i=1, 2
          pc(i+5)=pc(i+5)+phi(3,3,n)*pk(1+i)
        enddo
      enddo


      ENDSUBROUTINE

     !============================================================ 


