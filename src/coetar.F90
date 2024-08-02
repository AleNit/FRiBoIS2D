      
     !It computes the coefficients for the integration in directions
     !with non-uniform grid spacing, Y, Z
      
      SUBROUTINE coetar

      USE param
      
      implicit none
     !---------------------------------------------------------- 
      integer :: km,kc,kp
      integer :: jm,jc,jp
      integer :: inwn, inws
      real:: a33,a33m,a33p
      real:: a22,a22m,a22p
     !---------------------------------------------------------- 

      if (inslwn) then
        inwn=1    !no-slip condition top wall
      else 
        inwn=0       
      endif

      if (inslws) then
        inws=1    !no-slip condition bottom wall
      else
        inws=0
      endif

      
     !-------------------------------- Y direction 
      am2j(1)=0.0
      ap2j(1)=0.0
      ac2j(1)=1.0
      am2j(n2)=0.0
      ap2j(n2)=0.0
      ac2j(n2)=1.0
      
      do jc=2,n2m
       jm=jc-1
       jp=jc+1
       a22=dx2q/g2rc(jc)
       a22p=1.0/g2rm(jc)
       a22m=1.0/g2rm(jm)
       ap2j(jc)=a22*a22p
       am2j(jc)=a22*a22m
       ac2j(jc)=-(a22*a22p+a22*a22m)
      enddo


     !dens sidewall adiabatic
      do jc=2,n2m-1
       jp=jpv(jc)
       a22=dx2q/g2rm(jc)
       a22p= a22/g2rc(jp)
       a22m= a22/g2rc(jc)
       ap3j(jc)=a22p
       am3j(jc)=a22m
       ac3j(jc)=-(a22p+a22m)
       apscj(jc)=a22p
       amscj(jc)=a22m
       acscj(jc)=-(a22p+a22m)
      enddo

      jc=n2m
      jp=jc+1
      jm=jc-1
      a22=dx2q/g2rm(jc)
      a22p=+a22/g2rc(jp)
      a22m=+a22/g2rc(jc)
      apscj(jc)=0.0
      amscj(jc)=a22m
      acscj(jc)=-a22m
      am3j(jc)=a22m
      ap3j(jc)=2.0*inwn*a22p
      ac3j(jc)= -(a22m+inwn*a22p*2.0)

      jc=1
      jp=jc+1
      jm=jc-1
      a22=dx2q/g2rm(jc)
      a22p= +a22/g2rc(jp)
      a22m= +a22/g2rc(jc)
      apscj(jc)=a22p
      amscj(jc)=0.0
      acscj(jc)=-a22p
      ap3j(jc)=a22p
      am3j(jc)=2.0*a22m*inws
      ac3j(jc)=-(a22p+2.0*inws*a22m)

     
     !-------------------------------- Z direction 

     !coefficients for differentiation along Z
     !c means centered that is at k location
      am3ck(1)=0.0
      ap3ck(1)=0.0
      ac3ck(1)=1.0
      am3ck(n3)=0.0
      ap3ck(n3)=0.0
      ac3ck(n3)=1.0
      do kc=2, n3m
        km=kc-1
        kp=kc+1
        a33=dx3q/g3rc(kc)
        a33p=1.0/g3rm(kc)
        a33m=1.0/g3rm(km)
        ap3ck(kc)=a33*a33p
        am3ck(kc)=a33*a33m
        ac3ck(kc)=-(ap3ck(kc)+am3ck(kc))
      enddo

     !s means staggered that is at k+1/2 location 
      do kc=2,n3m-1
        kp=kc+1
        km=kc-1
        a33=dx3q/g3rm(kc)
        a33p= +a33/g3rc(kp)
        a33m= +a33/g3rc(kc)
        ap3sk(kc)=a33p
        am3sk(kc)=a33m
        ac3sk(kc)=-(ap3sk(kc)+am3sk(kc))
      enddo
    
      kc=1
      kp=kc+1
      a33=dx3q/g3rm(kc)
      a33p=+a33/g3rc(kp)
      a33m=+a33/g3rc(kc)
      ap3sk(kc)=a33p
      am3sk(kc)=0.0
      ac3sk(kc)=-(a33p+a33m*2.0)

      kc=n3m
      kp=kc+1
      a33=dx3q/g3rm(kc)
      a33p=+a33/g3rc(kp)
      a33m=+a33/g3rc(kc)
      am3sk(kc)=a33m
      ap3sk(kc)=0.0
      ac3sk(kc)=-(a33m+a33p*2.0)


      ENDSUBROUTINE

