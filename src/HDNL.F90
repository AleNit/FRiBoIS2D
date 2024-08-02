

     !compute convective terms 

      SUBROUTINE hdnl2
      
      USE param
      USE local_arrays, ONLY: q2,q3,dph
      USE mpi_param, ONLY: kstart,kend
      USE outflow_vars

      implicit none
     !----------------------------------------------------------
      integer :: kc,kp,jp,jm,jc,ic,im,ip
      integer :: kmm,kpp
      integer :: kstartp, kendp
      real    :: h22,h23,udx1,h21
     !----------------------------------------------------------


      if (kstart.eq.1) then
       kstartp=2
      else
       kstartp=kstart
      end if

      if (kend.eq.n3m) then
       kendp=n3m-1
      else
       kendp=kend
      end if


     !h term for the q2 momentum equation at i+1/2,j,k+1/2
      do kc=kstart, kend
        kmm=kmv(kc)
        kpp=kpv(kc)
        kp=kc+1

        do jc=2, n2m
          jm=jc-1
          jp=jc+1

           !h22=(q2*dq2/dy)/dt
            h22=((q2(jp,kc)+q2(jc,kc))      &     
          *(q2(jp,kc)+q2(jc,kc))             &
          -(q2(jm,kc)+q2(jc,kc))             &
          *(q2(jm,kc)+q2(jc,kc)))*udx2c(jc)*0.25

           dph(jc,kc)=-h22
        
        enddo
      enddo

      do kc=kstartp, kendp
        kmm=kmv(kc)
        kpp=kpv(kc)
        kp=kc+1
        do jc=2, n2m
          jm=jmv(jc)
            
         !h23=(q3*dq2/dz)/dt 
          h23=((q3(jc,kp)+q3(jm,kp))*      &
               (q2(jc,kpp)+q2(jc,kc))      &
              -(q3(jc,kc)+q3(jm,kc))*      &
              (q2(jc,kc)+q2(jc,kmm)))*udx3m(kc)*0.25

          dph(jc,kc)= dph(jc,kc)-h23

        enddo
      enddo


     !OUTFLOW
      if (kend.eq.n3m) then
        kc=n3m
        kmm=kc-1
        do jc=2, n2m
          jm=jmv(jc)
          
         !h23=(q3*dq2/dz)/dt 
          h23=(qb2n(jc)*(qb3n(jc)+qb3n(jm))   &
              -(q3(jc,kc)+q3(jm,kc))*            &
              (q2(jc,kc)+q2(jc,kmm))*0.5)*udx3m(kc)*0.5
          
          dph(jc,kc)=dph(jc,kc)-h23
      
        enddo
      endif

  
     !INFLOW
      if (kstart.eq.1) then
        kc=1
        kp=kc+1
        do jc=2, n2m
          jm=jmv(jc)
            
         !h23=(q3*dq2/dz)/dt 
          h23=((q3(jc,kp)+q3(jm,kp))*(q2(jc,kc)+    &
          q2(jc,kp))*0.5-qb2s(jc)*(qb3s(jc)+        &
          qb3s(jm)))*udx3m(kc)*0.5

          dph(jc,kc)=dph(jc,kc)-h23

        enddo    
      endif
      
      
      ENDSUBROUTINE

     !===================================================================

      SUBROUTINE hdnl3

      USE param
      USE local_arrays, ONLY: q2, q3, qcap
      USE mpi_param, ONLY: kstart, kend
      
      implicit none
     !----------------------------------------------------------
      integer :: jc, kc
      integer :: kmm, km, kp, jp, jmm, jpp, ic, imm, ipp
      real    :: h32,h33,h31
      real    :: udx1
      integer :: kstartp, kendp
     !----------------------------------------------------------


      if(kstart.eq.1) then
       kstartp=2
      else
       kstartp=kstart
      endif

      if (kend.eq.n3m) then
       kendp=n3m-1
      else
       kendp=kend
      end if


      do kc=kstartp, kendp
        kmm=kc-1
        km=kmv(kc)
        kp=kc+1

        do jc=1, n2m
          jmm=jmv(jc)
          jpp=jpv(jc)
          jp=jc+1

         !h32=(q2*dq3/dy)/dt
          h32=(((q2(jp,kc)+q2(jp,km))         &
              *(q3(jpp,kc)+q3(jc,kc)))        &
              -((q2(jc,kc)+q2(jc,km))         &
              *(q3(jc,kc)+q3(jmm,kc))))*udx2m(jc)*0.25

         !h33=(q3*dq3/dz)/dt
          h33=((q3(jc,kp)+q3(jc,kc))*         &
              (q3(jc,kp)+q3(jc,kc))           &
              -(q3(jc,kc)+q3(jc,kmm))*        &
              (q3(jc,kc)+q3(jc,kmm)))*udx3c(kc)*0.25

          qcap(jc,kc)=-(h32+h33)

        enddo
      enddo

      
      ENDSUBROUTINE

     !===================================================================

      SUBROUTINE hdnlro

      USE param
      USE local_arrays, ONLY: q2,q3,hro,dens
      USE mpi_param, ONLY: kstart,kend
      
      implicit none
     !----------------------------------------------------------
      integer :: jc,kc,ic
      integer :: kpp,km,kp,jp,jmm,jpp,ip,imm,ipp
      real    :: h32,h33,udx2,udx1,h31
     !----------------------------------------------------------


     !h term for the q3 momentum equation at i+1/2,j+1/2,k

!      udx1=dx1*0.5
      udx2=dx2*0.5

      do kc=kstart,kend
        km=kmv(kc)
        kpp=kpv(kc)
        kp=kc+1

        do jc=1,n2m
          jp=jc+1
          jmm=jmv(jc)
          jpp=jpv(jc)
           
         !d(rho q_r)/dt        
          h32=(q2(jp,kc)*(dens(jpp,kc)+dens(jc,kc))-   &
           q2(jc,kc)*(dens(jc,kc)+dens(jmm,kc)))*udx2

         !d(rho q_x)/dt
          h33=(q3(jc,kp)*(dens(jc,kpp)+dens(jc,kc))-   &
          q3(jc,kc)*(dens(jc,kc)+dens(jc,km)))*udx3m(kc)*0.5d0
          
          hro(jc,kc)=-(h32+h33)

        enddo
      enddo

      ENDSUBROUTINE

