
     !It computes the solenoidal velocity field at time-step
     !n+1, once the Poisson equation for pressure has been solved
      
      SUBROUTINE updvp
      
      USE param
      USE local_arrays,  ONLY: q2, q3, dph
      USE mpi_param,     ONLY: kstart, kend

      use mpih
      
      implicit none
     !----------------------------------------------------------
      integer :: jc, jm, kc, km, ic, im, kstartp
      real    :: usukm, usujm, udx1, locdph
     !----------------------------------------------------------

      if (kstart.eq.1) then
        kstartp=2
      else
        kstartp=kstart
      endif

     
     !Correct q2
      do kc=kstart, kend
        do jc=1, n2m
          jm=jmv(jc)
          usujm = al*dt*udx2c(jc)
          
          locdph=dph(jc,kc)
          q2(jc,kc)=q2(jc,kc)-(locdph-dph(jm,kc))*usujm
       
        enddo
      enddo


!     !re-set BCs... just for post-processing purposes, this leads to
!     !high divergence values for bodies in the wall vicinity
!
!      q2(:,1,:)=0.0
!      q2(:,n2,:)=0.0


     !Correct q3
      do kc=kstartp, kend
        km=kc-1
        usukm = al*dt*udx3c(kc)
        do jc=1,n2m
          
          q3(jc,kc)=q3(jc,kc)-(dph(jc,kc)-    &
                  dph(jc,km))*usukm
         
        enddo
      enddo

      
      ENDSUBROUTINE

     !===================================================================
     
     !It computes the pressure field at time-step n+1, once
     !the Poisson equation has been solved. The pressure
     !update formula is:
     ! p^(n+1) = p^(n) + phi^(n+1) - beta * Nabla^2 phi^(n+1)
      
      SUBROUTINE prcalc
      
      USE param
      USE local_arrays,  ONLY: pr, dph
      USE mpi_param,     ONLY: kstart, kend
      
      implicit none
     !----------------------------------------------------------
      integer :: kp, km, jm, jp, jc, kc, ic, ip, im
      real    :: be
      real :: ammk, acck, appk, ammj, accj, appj
     !----------------------------------------------------------

      be=al*beta

      do kc=kstart, kend
        kp=kpv(kc)
        km=kmv(kc)
        ammk=amphk(kc)
        acck=acphk(kc)
        appk=apphk(kc)

        do jc=1, n2m
          jm=jmv(jc)
          jp=jpv(jc)
          ammj=amphj(jc)
          accj=acphj(jc) 
          appj=apphj(jc)
          
          pr(jc,kc)=pr(jc,kc)+dph(jc,kc)-    &
            be*((dph(jp,kc)*appj+dph(jc,kc)*accj+dph(jm,kc)*ammj)+    &
                (dph(jc,kp)*appk+dph(jc,kc)*acck+dph(jc,km)*ammk))    

        enddo
      enddo


      ENDSUBROUTINE
     
     !===================================================================

     !subtract mean integral pressure from pr to get
     !zero mean value pressure field

      SUBROUTINE pr_submean

      USE param
      USE local_arrays,  ONLY: pr
      USE mpi_param,     ONLY: kstart, kend
      USE mpih

      implicit none
     !----------------------------------------------------------
      integer :: i, j, k, ip, jp, kp
      real :: mypravg, pravg
     !----------------------------------------------------------

      mypravg=0.0

      do k=kstart, kend
        kp=k+1
        do j=1, n2m
          jp=j+1
          mypravg=mypravg+pr(j,k)*(zz(kp)-zz(k))*(rc(jp)-rc(j)) 
        enddo
      enddo
            
      CALL MPI_ALLREDUCE(mypravg,pravg,1,MDP,MPI_SUM, &
               MPI_COMM_WORLD,ierr)

      pravg=pravg/(rext2*alx3)


      pr=pr-pravg


      ENDSUBROUTINE

 


