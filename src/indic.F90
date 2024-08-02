      
     !It computes local grid indices 
      
      SUBROUTINE indic
      
      USE param
      
      implicit none
     !----------------------------------------------------------
      integer :: jc, kc, ic
     !----------------------------------------------------------


     !---------------------- Z direction (3) - inflow/outflow
      do kc=1, n3m
        kmv(kc)=kc-1
        kpv(kc)=kc+1
        if (kc==1) kmv(kc)=kc
        if (kc==n3m) kpv(kc)=kc
      enddo
      do kc=1, n3m
        kpc(kc)=kpv(kc)-kc
        kmc(kc)=kc-kmv(kc)
        kup(kc)=1-kpc(kc)
        kum(kc)=1-kmc(kc)
      enddo


     !---------------------- Y direction (2) - wall
      do jc=1, n2
        jmv(jc)=jc-1
        jpv(jc)=jc+1
        if (jc.eq.1) jmv(jc)=jc
        if (jc.eq.n2m) jpv(jc)=jc
        if(jc.eq.n2) jpv(jc)=jc   
      enddo

      do jc=1, n2m
        jpc(jc)=jpv(jc)-jc
        jmc(jc)=jc-jmv(jc)
      enddo

      
      ENDSUBROUTINE

