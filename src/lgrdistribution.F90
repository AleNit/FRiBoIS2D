
     !find Eulerian cells containing each Lagrangian marker 

      SUBROUTINE lgrdistribution(ns)

      USE mpih
      USE utils_math
      USE param, ONLY: n2,n3,n2m,n3m,rext2,alx3,rm,rc,      &
            zm,zz,rext2,alx3
      USE mls_param
      USE rbm

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: ns
     !---------------------------------------------------------- 
      integer :: i, j, ind
      real :: yatmp(n2),zatmp(n3)
      real, dimension(2) :: lgrpos
      real :: rb_ymax,rb_zmax,rb_ymin,rb_zmin
     !----------------------------------------------------------

      
      if (myid==0) then 

       !check domain boundaries  
        rb_ymax=maxval(lmcc(:,1))
        rb_zmax=maxval(lmcc(:,2))
        rb_ymin=minval(lmcc(:,1))
        rb_zmin=minval(lmcc(:,2))

        if (rb_ymax>rext2.or.rb_ymin<0.0.or.    &
            rb_zmax>alx3.or.rb_zmin<0.0) then
          CALL stopcase
          stop
        endif 


        do i=1, nlm

          lgrpos=lmcc(i,:) 

         !forcing direction 1 
          yatmp=abs(rc-lgrpos(1))
          ind=minloc(yatmp,1)            
          lmind(i,1,1)=ind

          if (lgrpos(2)<zm(1).and.lgrpos(2)>0.0) then
            ind=1
          elseif (lgrpos(2)>zm(n3m).and.lgrpos(2)<alx3) then
            ind=n3m
          else
            zatmp=abs(zm-lgrpos(2))
            ind=minloc(zatmp,1)            
          endif
          lmind(i,1,2)=ind


         !forcing direction 2
          if (lgrpos(1)<rm(1).and.lgrpos(1)>0.0) then
            ind=1
          elseif (lgrpos(1)>rm(n2m).and.lgrpos(1)<rext2) then
            ind=n2m
          else
            yatmp=abs(rm-lgrpos(1))
            ind=minloc(yatmp,1)
          endif   
          lmind(i,2,1)=ind

          zatmp=abs(zz-lgrpos(2))
          ind=minloc(zatmp,1)
          lmind(i,2,2)=ind

        
        enddo   

      endif   
      
      CALL MPI_BCAST(lmind,nlm*4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
 

      ENDSUBROUTINE


