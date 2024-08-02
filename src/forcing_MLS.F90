
     !It forces boundary conditions on the fluid domain interpolating
     !the boundary values from a set of lagrangian markers
      
      SUBROUTINE forcing_MLS(substep)

      USE utils_math
      USE param,  ONLY: time, dimSD, dt, rc, zz, &
                        n2, n3, n2m, n3m,  &
                        nsst, rm, zm
      USE mls_param
      USE mpih
      USE mpi_param,  ONLY: kstart, kend
      USE rbm

      USE local_arrays, ONLY: q2f,q3f 

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: substep
     !---------------------------------------------------------- 
      integer :: i, j, k, idir, nlu, nlv, n, np, ktag
      integer :: ind, fs, kstartp, nlt_pa
      integer :: nz2, nz3
      real :: epu_l, epurn1, epurn2, epurn3
      real :: vol, lgr_h, dvel, areal
      integer, dimension(numtasks) :: rcounts_i, rcounts, displs
      integer, dimension(2) :: indse
      real, dimension(2) :: pos
      real, dimension(2,dimSD) :: SDc, SDd
      integer, dimension(2,2,dimSD) :: eulin
      real, dimension(:,:), allocatable :: epu, epur
     !----------------------------------------------------------
     

      do fs=1, forst  !---------- for iterative boundary conditions enforcement
      
      fey=0.0
      fez=0.0 

      if (substep==1.and.fs==1) then
        ktag=1
        ftag2=0
        ftag3=0         
        phipred=0.0      
        cl=0.0
        phi=0.0
      endif  
      
      
      do n=1, nlm   !---------- loop over lagrangian markers

       !check if the lagrangian marker falls within the process domain
        pos=lmcc(n,:)
        areal=lmA(n)

        if (zz(kstart)<pos(2).and.pos(2)<zz(kend+1)) then

          if (substep==1.and.fs==1) then

            do idir=1,2

              indse=lmind(n,idir,:) 
   
             !find indices and metrics of the support domain 
              CALL MLS_SD_ind(idir,indse,dimSD,SDc,SDd,eulin(idir,:,:))

             !compute transfer operator array phi 
              CALL MLS_SF(pos,indse,SDc,SDd,dimSD,phi(n,idir,:))
               
             !compute coefficient for momentum balance
              dvel=0.0
              lgr_h=0.0        
              do i=1, dimSD
                vol=SDd(1,i)*SDd(2,i)
                dvel=dvel+phi(n,idir,i)*vol
                lgr_h=lgr_h+1.0/2.0*phi(n,idir,i)*   &
                      (SDd(1,i)+SDd(2,i))
              enddo
    
              cl(n,idir)=areal*lgr_h/dvel              

            enddo

            
           !tag cells involved in the forcing procedure 
           !consider just one direction
            do i=1, dimSD
              ftag2(ktag)=eulin(1,1,i)
              ftag3(ktag)=eulin(1,2,i)
              ktag=ktag+1
            enddo

            eulins(n,:,:,:)=eulin

          endif

         
          CALL interp_spread(eulins(n,:,:,:),phi(n,:,:),lmvel(n,:),  &
              cl(n,:))

          endif  !---------- end marker selection 

      enddo   !---------- end loop over lagrangian markers

     
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if (substep==1.and.fs==1) then
        CALL MPI_ALLREDUCE(phi,phipred,   &
                 nlm*2*dimSD,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      endif

       
      CALL update_add_lower_ghost(n2,fey) 
      CALL update_add_lower_ghost(n2,fez) 
                                                
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)     
                                                
      CALL update_add_upper_ghost(n2,fey) 
      CALL update_add_upper_ghost(n2,fez) 

      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

      enddo   !--------- end loop for iterative bcs enforcement

      
      ENDSUBROUTINE
     
     !===================================================================

     !It performs field interpolation and spreading at the actual 
     !lagrangian marker
          
      SUBROUTINE interp_spread(eul,phi,vel,cl)
      
      USE param,  ONLY: dt, dimSD, n2m, n3m 
      USE local_arrays
      USE mls_param,  ONLY: fey, fez
      USE mpih

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: eul(2,2,dimSD)
      real, intent(in) :: phi(2,dimSD), vel(2), cl(2)
     !----------------------------------------------------------
      integer :: n
      real :: ucy, ucz, fcy, fcz
      integer :: tmp(2)
      real, dimension(dimSD) :: uk
     !----------------------------------------------------------
     

      ucy=0.0
      do n=1, dimSD
        uk(n)=q2f(eul(1,1,n),eul(1,2,n))
        ucy=ucy+phi(1,n)*uk(n)
      enddo
      fcy=(vel(1)-ucy)/dt

      do n=1, dimSD
        tmp=[eul(1,1,n),eul(1,2,n)] 
        fey(tmp(1),tmp(2))=fey(tmp(1),tmp(2))+  &
            cl(1)*phi(1,n)*fcy
      enddo


      ucz=0.0
      do n=1, dimSD
        uk(n)=q3f(eul(2,1,n),eul(2,2,n))
        ucz=ucz+phi(2,n)*uk(n)
      enddo
      fcz=(vel(2)-ucz)/dt

      do n=1, dimSD
        tmp=[eul(2,1,n),eul(2,2,n)]
        fez(tmp(1),tmp(2))=fez(tmp(1),tmp(2))+  &
               cl(2)*phi(2,n)*fcz
      enddo


      ENDSUBROUTINE

     !===================================================================

