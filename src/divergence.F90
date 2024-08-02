
      MODULE divergence

      USE param
      USE local_arrays,  ONLY: q2,q3,dph
      USE mpi_param,     ONLY: kstart,kend

      CONTAINS
     
     !===================================================================

     !It computes the local divergence of the non-solenoidal velocity 
     !field (qcap) and the pressure correction term (dph)
      
      SUBROUTINE divg

      USE mpih

      implicit none
     !----------------------------------------------------------
      integer :: jc, jp, kc, kp, ic, ip
      real    :: usdtal, dqcap   
     !----------------------------------------------------------

     !coefficient of div(u*) in the Poisson eqn  
      usdtal = 1.0/(dt*al)

      do kc=kstart, kend
        kp=kc+1

        do jc=1, n2m
          jp=jc+1            
          
          dqcap=(q2(jp,kc)-q2(jc,kc))*udx2m(jc)   &
               +(q3(jc,kp)-q3(jc,kc))*udx3m(kc)   
          dph(jc,kc)=dqcap*usdtal
           
         enddo
      enddo

    
      ENDSUBROUTINE

     !===================================================================

     !It find the maximum value of flow divergence on the 
     !computational domain  
      
      SUBROUTINE divgck(qmax)
      
      USE mpih
      
      implicit none
     !----------------------------------------------------------
      real,intent(out) :: qmax
      integer :: jc,kc,kp,jp,ic,ip
      real    :: dqcap, my_qmax
     !----------------------------------------------------------
       
      my_qmax =-100.0
      qmax=0.0                                                     
      
      do kc=kstart, kend
        kp=kc+1

        do jc=1,n2m
          jp=jc+1

           dqcap=(q2(jp,kc)-q2(jc,kc))*udx2m(jc)+   &
                 (q3(jc,kp)-q3(jc,kc))*udx3m(kc)
           my_qmax = max(abs(dqcap),my_qmax)      
 
        enddo
      enddo

      CALL MPI_ALLREDUCE(my_qmax,qmax,1,MDP,MPI_MAX,  &
                         MPI_COMM_WORLD,ierr)
      
      ENDSUBROUTINE        

     !===================================================================

     !It find the location(s) of excessive divergence 
     !i.e. i,j,k where dqcap(i,j,k) > resid
      
      SUBROUTINE divgloc
      
      USE mpih
      
      implicit none
     !----------------------------------------------------------
      integer :: jc, kc, kp, jp, ic, ip
      integer :: pldivg(kend-kstart+1)
      real    :: dqcap
     !----------------------------------------------------------
      
      pldivg=0
      
      do kc=kstart, kend
        kp=kc+1
        pldivg(kc-kstart+1)=0

        do jc=1, n2m
          jp=jc+1
          
          dqcap= (q2(jp,kc)-q2(jc,kc))*udx2m(jc)    &
                +(q3(jc,kp)-q3(jc,kc))*udx3m(kc)
          
           if (abs(dqcap)>resid) then
             pldivg(kc-kstart+1)=1
           endif
        
        enddo

      enddo

      
      ENDSUBROUTINE         

     !===================================================================
     
      ENDMODULE
