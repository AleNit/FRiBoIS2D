     
     !assemble right-hand side of the momentum equation

      SUBROUTINE invtr2

      USE mpih
      USE param
      USE local_arrays,  ONLY: q2, pr, rhs2, dph, ru2, q2f
      USE mpi_param, ONLY: kstart,kend
      USE outflow_vars

      implicit none
     !----------------------------------------------------------
      integer :: jc,kc,km,kp,jp,jm,ic,im,ip
      integer :: n,i,j,k,ie,je,ke
      integer :: kstartp, kendp
      real    :: amm,app,acc
      real    :: ammj,appj,accj
      real    :: dcq2,dpx22,q2e
      real    :: d22q2,d33q2,d11q2
      real    :: alre,udx1q,udx2
      real    :: checksum
      real    :: codea, codeb
     !----------------------------------------------------------
     

      alre=al/ren

     !compute the rhs of the factored equation
     !everything at i+1/2,j,k+1/2

     !points inside the flowfield
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

      do kc=kstartp, kendp
        km=kmv(kc)
        kp=kpv(kc)
        amm=am3sk(kc)
        acc=ac3sk(kc)
        app=ap3sk(kc)

        do jc=2, n2m
          jm=jc-1
          jp=jc+1
          ammj=am2j(jc)
          accj=ac2j(jc)
          appj=ap2j(jc)
          udx2=al*udx2c(jc)
          
          d22q2=q2(jp,kc)*appj+q2(jc,kc)*accj    &
                +q2(jm,kc)*ammj
          
          d33q2=q2(jc,kp)*app+q2(jc,kc)*acc      &
                +q2(jc,km)*amm

          dcq2=d22q2+d33q2

         !pressure gradient
          dpx22=(pr(jc,kc)-pr(jm,kc))*udx2

         !collect viscous, pressure, and convective terms and assemble RHS
          rhs2(jc,kc)=(ga*dph(jc,kc)+ro*ru2(jc,kc)   &
                        +alre*dcq2-dpx22)*dt

          ru2(jc,kc)=dph(jc,kc)         
          
        enddo
      enddo


      if (kstart.eq.1) then
        kc=1
        kp=kc+1
        codea = dx3q/(g3rc(kp)*g3rm(kc))
        codeb = dx3q/(g3rc(kc)*g3rm(kc))
        do jc=2, n2m
          jp=jc+1
          jm=jc-1
          ammj=am2j(jc)
          accj=ac2j(jc)
          appj=ap2j(jc)
          
          d22q2=q2(jp,kc)*appj+q2(jc,kc)*accj    &
                 +q2(jm,kc)*ammj

          d33q2=((q2(jc,kp)-q2(jc,kc))*codea-    & 
                (q2(jc,kc)-qb2s(jc) )*codeb*2.0)

          dcq2=d22q2+d33q2

          dpx22=(pr(jc,kc)-pr(jm,kc))*udx2c(jc)*al

          rhs2(jc,kc)=+(ga*dph(jc,kc)+ro*ru2(jc,kc)     &
                        +alre*dcq2-dpx22)*dt

          ru2(jc,kc)=dph(jc,kc)
         
        enddo
      endif


      if (kend.eq.n3m) then
        kc=n3m
        km=kc-1
        kp=kc+1
        codea = dx3q/(g3rc(kp)*g3rm(kc))
        codeb = dx3q/(g3rc(kc)*g3rm(kc))
        do jc=2, n2m
          jm=jc-1
          jp=jc+1
          ammj=am2j(jc)
          accj=ac2j(jc)
          appj=ap2j(jc)
          
          d22q2=q2(jp,kc)*appj+q2(jc,kc)*accj    &
                +q2(jm,kc)*ammj

          d33q2=((qb2n(jc)-q2(jc,kc))*codea*2.0-    &
                (q2(jc,kc)-q2(jc,km))*codeb ) 

          dcq2=d22q2+d33q2

          dpx22=(pr(jc,kc)-pr(jm,kc))*udx2c(jc)*al

          rhs2(jc,kc)=+(ga*dph(jc,kc)+ro*ru2(jc,kc)     &
                        +alre*dcq2-dpx22)*dt

          ru2(jc,kc)=dph(jc,kc)
        
        enddo
      endif
     
     
     !compute preliminary velocity field
      do kc=kstart, kend
        do jc=2, n2m
          q2f(jc,kc)=q2(jc,kc)+rhs2(jc,kc)
        enddo
      enddo


      ENDSUBROUTINE

     !===================================================================

      SUBROUTINE invtr3
      
      USE param
      USE local_arrays, ONLY: q3, qcap, pr, ru3, rhs3, q3f
      USE mpi_param, ONLY: kstart, kend
      USE outflow_vars
      
      USE mpih

      implicit none
     !----------------------------------------------------------
      integer :: jc,kc,km,kp,jp,jm,ic,ip,im
      integer :: i,j,k,ie,je,ke,n
      real    :: udx3
      real    :: dq32,dq33,dcq3,dpx33,dq31
      real    :: app,acc,amm,q3e
      real    :: appj,accj,ammj
      real    :: alre,udx1q, checksum
      integer :: kstartp
     !----------------------------------------------------------

      if(kstart.eq.1) then
       kstartp=2
      else
       kstartp=kstart
      endif

      alre=al/ren

     
     !compute the rhs of the factored equation
     !everything at i+1/2,j+1/2,k
      do kc=kstartp, kend
        km=kc-1
        kp=kc+1
        udx3=al*udx3c(kc)
        amm=am3ck(kc)
        acc=ac3ck(kc)
        app=ap3ck(kc)
        
       !Moving lower wall
        jc=1
        jp=jc+1
        ammj=am3j(jc)
        accj=ac3j(jc)
        appj=ap3j(jc)
        
        dq32=q3(jp,kc)*appj+q3(jc,kc)*accj    &
             +q3lw*ammj

        dq33=q3(jc,kp)*app +q3(jc,kc)*acc     &
             +q3(jc,km)*amm

        dcq3=dq32+dq33

       !pressure gradient
        dpx33=(pr(jc,kc)-pr(jc,km))*udx3

       !collect for RHS assembly 
        rhs3(jc,kc)=(ga*qcap(jc,kc)+ro*ru3(jc,kc)  &
                       +alre*dcq3-dpx33)*dt 

        ru3(jc,kc)=qcap(jc,kc)
       
       !Inner points
        do jc=2, n2m-1
          jm=jmv(jc)
          jp=jpv(jc)
          ammj=am3j(jc)
          accj=ac3j(jc)
          appj=ap3j(jc)
          
          dq32=q3(jp,kc)*appj+q3(jc,kc)*accj     &
              +q3(jm,kc)*ammj

          dq33=q3(jc,kp)*app +q3(jc,kc)*acc      &
               +q3(jc,km)*amm

          dcq3=dq32+dq33

          dpx33=(pr(jc,kc)-pr(jc,km))*udx3

          rhs3(jc,kc)=(ga*qcap(jc,kc)+ro*ru3(jc,kc)     &
                       +alre*dcq3-dpx33)*dt 

          ru3(jc,kc)=qcap(jc,kc)
        
        enddo
       
     
       !Upper wall
        jc=n2m
        jm=jc-1
        ammj=am3j(jc)
        accj=ac3j(jc)
        appj=ap3j(jc)
        
        dq32=q3uw*appj+q3(jc,kc)*accj      &
            +q3(jm,kc)*ammj

        dq33=q3(jc,kp)*app +q3(jc,kc)*acc    &
                +q3(jc,km)*amm

        dcq3=dq32+dq33

        dpx33=(pr(jc,kc)-pr(jc,km))*udx3

        
        rhs3(jc,kc)=(ga*qcap(jc,kc)+ro*ru3(jc,kc)   &
                         +alre*dcq3-dpx33)*dt 

        ru3(jc,kc)=qcap(jc,kc)
     
      enddo

     
     !compute preliminary velocity field
      do kc=kstartp, kend
        do jc=1, n2m
          q3f(jc,kc)=q3(jc,kc)+rhs3(jc,kc)
        enddo
      enddo

      
      ENDSUBROUTINE
     
     !===================================================================


