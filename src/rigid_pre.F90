   
     !pre-processing operations for the rigid body dynamics

      SUBROUTINE rigid_pre(i_area)
     
      USE utils_math
      USE param,  ONLY: dt,dt_o,time,outfol,tins,irest
      USE rbm
      USE mpih
      USE mls_param, ONLY: cbody,ddmin

      implicit none
     !----------------------------------------------------------
      real, intent(out) :: i_area
     !----------------------------------------------------------
      integer :: i, j, k, l, z, c
      integer :: ll, lln
      integer :: load_dir, ndof_tot_old
      integer :: config, info
      integer :: tmpreadi, msl
      real :: f_area, time0
      real :: small
      integer, allocatable :: double(:), ms(:,:)
      real, allocatable :: rs(:), RR(:,:)
      character(50) :: newdir
      character(100) :: filename
      character(20) :: patch_ID

      real :: radius(2),posP(2),velP(2),accP(2),rho(2)
      real :: crprod1(2),crprod2(2),crprod3(2)
      real :: rmat(2,2), rmati(2,2)
     !---------------------------------------------------------- 

      i_area=0.0
      
      if (myid==0) then
      

      CALL cpu_time(tins(1))

         
       CALL read_input_par_RBM(rbden,rba,    &
             cmposi,cmroti,rbgravF,rbmcon,cbody,kk,    &
             cc,kk3,rtcoup,rrr,ipert)

       
      !initialize center of mass kinematics 
       cmp_n(1:2)=cmposi
       cmp_n(3)=cmroti
       cmv_n=0.0
       cma_n=0.0


       kinori(1,1:2)=cmposi
       kinori(1,3)=cmroti
       kinori(2,:)=cmv_n
       kinori(3,:)=cma_n

      
      !overwrite initial perturbation
       cmv_n=ipert
       if (rtcoup)  cmv_n(3)=cmv_n(1)/rrr       !kinematic coupling (rotation-translation)


      !over-write rigid body kinematics in case of restart
       if (irest) then
   
         filename='./input_FSI/restart/master_rigid.dat'
         open(447,file=trim(filename),status='old')
         read(447,*) time
         read(447,*) dt,dt_o
         read(447,*) cmp_n
         read(447,*) cmv_n
         read(447,*) cma_n
         close(447)
   
       endif


       cmp_np1=cmp_n
       cmv_np1=cmv_n
       cma_np1=cma_n


      !initialize Lagrangian markers kinematics
       CALL get_rotmat(cmp_n(3),rmat,rmati)

       do i=1, nlm
            
          rho=lmcc(i,:)-cmp_n(1:2)
   
          posP=cmp_np1(1:2)+matmul(rmat,rho) 
   
          crprod1(1)=-cmv_np1(3)*rho(2)
          crprod1(2)=cmv_np1(3)*rho(1)
          velP=cmv_np1(1:2)+crprod1
   
          crprod2(1)=-cma_np1(3)*rho(2)
          crprod2(2)=cma_np1(3)*rho(1)
          crprod3(1)=-cmv_np1(3)*crprod1(2)
          crprod3(2)=cmv_np1(3)*crprod1(1)
          accP=cma_np1(1:2)+crprod2+crprod3
   
          lmcc(i,:)=posP
          lmvel(i,:)=velP
          lmacc(i,:)=accP
   
       enddo


      !apply initial rotation to vertices and recompute normal arrays
       do i=1, nlm+1
         radius=vert(i,:)-cmposi 
         posP=cmp_np1(1:2)+matmul(rmat,radius)
         vert(i,:)=posP
       enddo   

      CALL createnorm(nlm,vert,-1.0,lmnor)


      CALL cpu_time(tins(2))

      endif    !end scalar section

      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      
      CALL MPI_BCAST(cbody,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ddmin,1,MDP,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(time,1,MDP,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dt,1,MDP,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dt_o,1,MDP,0,MPI_COMM_WORLD,ierr)

      CALL MPI_BCAST(vert,(nlm+1)*2,MDP,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lmcc,nlm*2,MDP,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lmvel,nlm*2,MDP,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lmacc,nlm*2,MDP,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lmnor,nlm*2,MDP,0,MPI_COMM_WORLD,ierr)


      ENDSUBROUTINE


