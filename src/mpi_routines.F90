
     !get start and end indices for each task
      
      SUBROUTINE block(n, p, irank, istart, iend, blcsz)
      
      implicit none
     !----------------------------------------------------------
      integer,intent(in) :: n,p,irank
      integer,intent(out) :: istart,iend
      integer :: i
      integer,dimension(0:p-1),intent(out) :: blcsz
     !----------------------------------------------------------
     
      do i=0,p-1
        blcsz(i) = floor(real((n+p-i-1)/p))
      enddo
      
      istart = sum(blcsz(0:irank))-blcsz(irank)+1
      iend = istart+blcsz(irank)-1

      ENDSUBROUTINE
      
     !===================================================================
     
      SUBROUTINE mpi_workdistribution
      
      USE param
      USE mpih 
      USE mpi_param
      
      implicit none
     !----------------------------------------------------------
      integer :: i, n1
     !----------------------------------------------------------

     
     !define number of ghost cells from MLS support domain size
      lvlhalo=3

      if(.not. allocated(counti)) allocate(counti(0:numtasks-1))
      if(.not. allocated(countj)) allocate(countj(0:numtasks-1))
      if(.not. allocated(countjp)) allocate(countjp(0:numtasks-1))
      if(.not. allocated(countk)) allocate(countk(0:numtasks-1))
   
     !For PERIODIC pressure solver
      CALL block(n2+1, numtasks, myid, jstartp, jendp, countjp)
      djp=jendp-jstartp+1
    
!      CALL block(n1+1, numtasks, myid, istart, iend, counti)
!      di=iend-istart+1
    
      CALL block(n2m, numtasks, myid, jstart, jend, countj)
      dj=jend-jstart+1
     
      CALL block(n3m, numtasks, myid, kstart, kend, countk)
      dk=kend-kstart+1


#ifdef DEBUG
      if (myid==0) then
        print*,''
        print*, 'Indices for MPI job distribution'
      endif
      do i=0, numtasks-1
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (myid==i) then
          print*,''
          print*, ' ------------------------------ Task ID:', myid
          write(*,*) " jstart:   ", jstart
          write(*,*) " jend:     ", jend
          write(*,*) " jstartp:  ", jstart
          write(*,*) " jendp:    ", jend
          write(*,*) " kstart:   ", kstart
          write(*,*) " kend:     ", kend
        endif
      enddo
#endif

      if (mod(n3m,numtasks)/=0) then
        if (myid==0) then
        print*,''
        print*, '---------------------------------'
        print*, ' The number of cells in direction z is not a '// &
                'multiple of the number of tasks...'
        print*, '---------------------------------'
        print*,''
        stop
        endif
      endif
   
      if (dj.lt.1) then
        write(6,*) 'process ',myid,' has work load <1 cell in j dir'
        write(6,*) 'Check grid dimensions and number of processes'
        CALL MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif
      if (dk.lt.1) then
        write(6,*) 'process ',myid,' has work load <1 cell in k dir'
        write(6,*) 'Check grid dimensions and number of processes'
        CALL MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif
    
      if (.not. allocated(offseti)) allocate(offseti(0:numtasks-1))
      if (.not. allocated(offsetjp)) allocate(offsetjp(0:numtasks-1))
      if (.not. allocated(offsetj)) allocate(offsetj(0:numtasks-1))
      if (.not. allocated(offsetk)) allocate(offsetk(0:numtasks-1))
                                             
!      offseti(:)=0
      offsetjp(:)=0
      offsetj(:)=0
      offsetk(:)=0
      do i=1,numtasks-1
!        offseti(i)= offseti(i-1) + counti(i-1)
        offsetjp(i)= offsetjp(i-1) + countjp(i-1)
        offsetj(i)= offsetj(i-1) + countj(i-1)
        offsetk(i)= offsetk(i-1) + countk(i-1)
      enddo
     
     
     !-------For MPI-IO--------------------------------
      mydata= n2*dk
      mydatam = n2m*dk

      n1=1
     
      if (myid.eq.numtasks-1)  mydata=n2*(dk+1)*n1
     
      if (.not. allocated(countf)) allocate(countf(0:numtasks-1))
      if (.not. allocated(offsetf)) allocate(offsetf(0:numtasks-1))
     
      CALL MPI_ALLGATHER(mydata,1,MPI_INTEGER,countf,1,MPI_INTEGER,  &
           MPI_COMM_WORLD,ierr)
    
      offsetf(:)=0
      do i=1,numtasks-1
        offsetf(i)= offsetf(i-1) + countf(i-1)
      enddo       
     !------------------------------------------------
      
      ENDSUBROUTINE

     !===================================================================

      SUBROUTINE update_both_ghosts(n2,q1,ks,ke)
      
      USE mpih
      
      implicit none
     !----------------------------------------------------------      
      integer, intent(in) :: ks,ke
      real,intent(inout) :: q1(n2,ks-lvlhalo:ke+lvlhalo)
      integer,intent(in) :: n2
      integer :: mydata
      integer :: my_down, my_up,tag
     !----------------------------------------------------------

      mydata= n2*lvlhalo

      my_down=myid-1
      
      my_up=myid+1

      if(myid .eq. 0) my_down=MPI_PROC_NULL
      if(myid .eq. numtasks-1) my_up=MPI_PROC_NULL

      tag=1
      CALL MPI_ISEND(q1(1,ke-lvlhalo+1), mydata, MDP, &
         my_up,tag,MPI_COMM_WORLD,req(1),ierr)
      
      CALL MPI_ISEND(q1(1,ks), mydata,  MDP, &
         my_down,tag,MPI_COMM_WORLD,req(2), ierr)
     
      CALL MPI_IRECV(q1(1,ks-lvlhalo), mydata,  MDP, &
         my_down,tag,MPI_COMM_WORLD,req(3),ierr)
     
      CALL MPI_IRECV(q1(1,ke+1), mydata,  MDP, &
         my_up, tag,MPI_COMM_WORLD,req(4),ierr)
     
      CALL MPI_Waitall(4,req,status,ierr)

      ENDSUBROUTINE

     !===================================================================

      SUBROUTINE update_upper_ghost_nohalo(n2,q1)
      
      USE mpih
      USE mpi_param, ONLY: kstart,kend,dk
      
      implicit none
     !---------------------------------------------------------- 
      integer,intent(inout) :: q1(n2,kstart-1:kend)
      integer,intent(in) :: n2
      integer :: mydata
      integer :: my_down, my_up,tag
     !---------------------------------------------------------- 
       
      mydata= n2
      
      my_down= myid-1
      
      my_up= myid+1

      if(myid .eq. 0) my_down= MPI_PROC_NULL
      if(myid .eq. numtasks-1) my_up= MPI_PROC_NULL
     
      tag=1
      
      CALL MPI_ISEND(q1(1,kstart-1), mydata, MDP, &
          my_down, tag, MPI_COMM_WORLD, req(1), ierr)
      
      CALL MPI_IRECV(q1(1,kend), mydata, MDP, &
          my_up,tag, MPI_COMM_WORLD, req(2), ierr)
      
      CALL MPI_Waitall(2,req,status,ierr)

      ENDSUBROUTINE

     !===================================================================

      SUBROUTINE update_upper_ghost(n2,q1)
      
      USE mpih
      USE mpi_param, only: kstart,kend,dk
      
      implicit none
     !---------------------------------------------------------- 
      real,intent(inout) :: q1(n2,kstart-lvlhalo:kend+lvlhalo)
      integer,intent(in) :: n2
      integer :: mydata
      integer :: my_down, my_up,tag
     !---------------------------------------------------------- 
       
      mydata= n2*lvlhalo
      
      my_down= myid-1
      
      my_up= myid+1

      if(myid .eq. 0) my_down= MPI_PROC_NULL
      if(myid .eq. numtasks-1) my_up= MPI_PROC_NULL
     
      tag=1
      
      CALL MPI_ISEND(q1(1,kstart), mydata, MDP, &
           my_down, tag, MPI_COMM_WORLD, req(1), ierr)
      
      CALL MPI_IRECV(q1(1,kend+1), mydata, MDP, &
           my_up,tag, MPI_COMM_WORLD, req(2), ierr)
       
      CALL MPI_Waitall(2,req,status,ierr)
    
      ENDSUBROUTINE

     !===================================================================

      SUBROUTINE update_lower_ghost(n2,q1)
      
      USE mpih
      
      USE mpi_param, only: kstart,kend,dk
      
      implicit none
     !---------------------------------------------------------- 
      real,intent(inout) :: q1(n2,kstart-lvlhalo:kend+lvlhalo)
      integer,intent(in) :: n2
      integer :: mydata
      integer :: my_down, my_up,tag
      integer :: ic,jc,kc
      real :: cksum
     !---------------------------------------------------------- 
       
      mydata=n2*lvlhalo
      
      my_down= myid-1
      
      my_up= myid+1

      if(myid .eq. 0) my_down= MPI_PROC_NULL
      if(myid .eq. numtasks-1) my_up= MPI_PROC_NULL
      
      tag=1
      
      CALL MPI_ISEND(q1(1,kend-lvlhalo+1), mydata,  MDP, &
           my_up, tag, MPI_COMM_WORLD, req(1), ierr)
      
      CALL MPI_IRECV(q1(1,kstart-lvlhalo), mydata,  MDP, &
           my_down,tag, MPI_COMM_WORLD, req(2), ierr)
       
      CALL MPI_Waitall(2,req,status,ierr)

      ENDSUBROUTINE

     !===================================================================

      SUBROUTINE update_add_upper_ghost(n2,q)
      
      USE mpih
      USE mpi_param, ONLY: kstart,kend,dk
      
      implicit none
     !---------------------------------------------------------- 
      real,intent(inout) :: q(n2,kstart-lvlhalo:kend+lvlhalo-1)
      real :: buf(n2,lvlhalo)
      integer,intent(in) :: n2
      integer :: mydata
      integer :: my_down, my_up,tag
      integer :: ic,jc
      real :: cksum
     !----------------------------------------------------------
       
      mydata= n2*lvlhalo
      
      my_down= myid-1
      
      my_up= myid+1

      buf=0.0d0

      if(myid .eq. 0) my_down= MPI_PROC_NULL
      if(myid .eq. numtasks-1) my_up= MPI_PROC_NULL
     
      tag=1

      CALL MPI_ISEND(q(1,kstart-lvlhalo),mydata,MDP, &
           my_down, tag, MPI_COMM_WORLD, req(1), ierr)
      
      CALL MPI_IRECV(buf(1,1), mydata, MDP, &
           my_up,tag, MPI_COMM_WORLD, req(2), ierr)

      CALL MPI_Waitall(2,req,status,ierr)

      do ic=1,lvlhalo 
        jc=kend-lvlhalo+ic
        if (jc.eq.kend) then
          q(:,jc) = buf(:,ic)
        else
          q(:,jc) = q(:,jc) + buf(:,ic)
        endif
      enddo

      ENDSUBROUTINE

     !=================================================================== 
      
      SUBROUTINE update_add_lower_ghost(n2,q)
      
      USE mpih
      USE mpi_param, only: kstart,kend,dk
      
      implicit none
     !---------------------------------------------------------- 
      real,intent(inout) :: q(n2,kstart-lvlhalo:kend+lvlhalo-1)
      real :: buf(n2,lvlhalo)
      integer,intent(in) :: n2
      integer :: mydata
      integer :: my_down, my_up,tag
      integer :: ic,jc
     !---------------------------------------------------------- 
       

      mydata= n2*lvlhalo
      
      my_down= myid-1
      
      my_up= myid+1

      buf=0.0d0

      if(myid .eq. 0) my_down= MPI_PROC_NULL
      if(myid .eq. numtasks-1) my_up= MPI_PROC_NULL
      
      tag=1
      
      CALL MPI_ISEND(q(1,kend),mydata,MDP, &
           my_up, tag, MPI_COMM_WORLD, req(1), ierr)

      CALL MPI_IRECV(buf(1,1), mydata, MDP, &
           my_down,tag, MPI_COMM_WORLD, req(2), ierr)

      CALL MPI_Waitall(2,req,status,ierr)

      do ic=1,lvlhalo
        jc=kstart+ic-2
        q(:,jc) = q(:,jc) + buf(:,ic)
      enddo

      ENDSUBROUTINE

     !===================================================================

      SUBROUTINE mpi_globalsum_double_arr(var,nvar)
      
      USE mpih
      
      implicit none
     !---------------------------------------------------------- 
      real,intent(inout),dimension(nvar) :: var
      real,dimension(nvar) :: var2
      integer,intent(in) :: nvar
     !---------------------------------------------------------- 

      CALL MPI_ALLREDUCE(var,var2,nvar,MPI_DOUBLE_PRECISION, &
            MPI_SUM,MPI_COMM_WORLD,ierr)          
      
      var = var2

      ENDSUBROUTINE

     !===================================================================

      SUBROUTINE mem_alloc
      
      USE mpih
      USE param, only: n2
      USE mpi_param
      USE local_arrays
     
      implicit none
     !---------------------------------------------------------- 
      integer :: merr(100), k
     !----------------------------------------------------------
 
      merr=0
      k=1
   
   
     !Arrays with ghost cells
      allocate(q2(1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr(k))
      k=k+1
      allocate(q3(1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr(k))
      k=k+1
      allocate(pr(1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr(k))
      k=k+1
      allocate(dens(1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr(k))
      k=k+1
      allocate(dph(1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr(k))
      k=k+1
       
      allocate(q2f(1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr(k))
      k=k+1
      allocate(q3f(1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr(k))
      k=k+1

     !storage for strong coupling
      allocate(q2_o(1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr(k))
      k=k+1
      allocate(q3_o(1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr(k))
      k=k+1
      allocate(pr_o(1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr(k))
      k=k+1
      allocate(dens_o(1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr(k))
      k=k+1

 

     !Arrays without ghost cells
      allocate(rhs(1:n2,kstart:kend), stat=merr(k))
      k=k+1
      allocate(ru2(1:n2,kstart:kend), stat=merr(k))
      k=k+1
      allocate(ru3(1:n2,kstart:kend), stat=merr(k))
      k=k+1
      allocate(qcap(1:n2,kstart:kend), stat=merr(k))
      k=k+1
      allocate(hro(1:n2,kstart:kend), stat=merr(k))
      k=k+1
      allocate(ruro(1:n2,kstart:kend), stat=merr(k))
      k=k+1
      allocate(rhs2(1:n2,kstart:kend), stat=merr(k))
      k=k+1
      allocate(rhs3(1:n2,kstart:kend), stat=merr(k))
 
     !storage for strong coupling1
      allocate(ru2_o(1:n2,kstart:kend), stat=merr(k))
      k=k+1
      allocate(ru3_o(1:n2,kstart:kend), stat=merr(k))


     !check successful allocation 
      if (any(merr/=0)) then
        write(6,*)"process  ",myid," failed to allocate memory"
      endif

          
      ENDSUBROUTINE

     !===================================================================
      
      SUBROUTINE mem_dealloc
      
      USE local_arrays
      USE mpi_param
      
      implicit none
      
!      if(allocated(q1)) deallocate(q1)
      if(allocated(q2)) deallocate(q2)
      if(allocated(q3)) deallocate(q3)
      
      if(allocated(qcap)) deallocate(qcap)
      
      if(allocated(dens)) deallocate(dens)
      if(allocated(pr)) deallocate(pr)
      if(allocated(hro)) deallocate(hro)
      
      if(allocated(rhs)) deallocate(rhs)
      
      if(allocated(dph)) deallocate(dph)
      
!      if(allocated(ru1)) deallocate(ru1)
      if(allocated(ru2)) deallocate(ru2)
      if(allocated(ru3)) deallocate(ru3)
      
      if(allocated(ruro)) deallocate(ruro)
      
      !---------------------------------------
      if(allocated(countj)) deallocate(countj)
      if(allocated(countk)) deallocate(countk)

      if(allocated(offsetj)) deallocate(offsetj)
      if(allocated(offsetk)) deallocate(offsetk)
      
      if(allocated(countf)) deallocate(countf)
      
      if(allocated(offsetf)) deallocate(offsetf)
      
   
      ENDSUBROUTINE

     !===================================================================

      SUBROUTINE mpi_write_continua(ns)
      
      USE param
      USE mpih
      USE mpi_param, ONLY: kstart,kend
      USE local_arrays, ONLY: q2, q3, pr, dens, ru2, ru3
      USE hdf5
      
      implicit none
     !----------------------------------------------------------
      integer :: ns
      integer :: hdf_error, comm, info, ndims
      integer :: i
      integer(HID_T) :: file_id, filespace, slabspace, memspace
      integer(HID_T) :: dset_q1, dset_q2, dset_q3, dset_dens, dset_pr
      integer(HID_T) :: plist_id, dset_grid, dspace_grid
      integer(HSIZE_T) :: dims(2), dims_grid(1)
      integer(HSIZE_T), dimension(2) :: data_count  
      integer(HSSIZE_T), dimension(2) :: data_offset 
      character(50) :: filnam1, filnam2, filnam3, filnam4, filnam5
      character(50) :: filnam6, filnam7, filnam8
      character(50) :: filnamgrid
     !----------------------------------------------------------

     !Sort out MPI definitions
      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

     !Form the name of the file
      filnam1 = trim(outfol)//'/restart/rest_dens.h5'
      filnam3 = trim(outfol)//'/restart/rest_q2.h5'
      filnam4 = trim(outfol)//'/restart/rest_q3.h5'
      filnam5 = trim(outfol)//'/restart/rest_pr.h5'
      filnam7 = trim(outfol)//'/restart/rest_ru2.h5'
      filnam8 = trim(outfol)//'/restart/rest_ru3.h5'

     !Set offsets and element counts
      ndims = 2

      dims(1)=n2
      dims(2)=n3m

      CALL h5screate_simple_f(ndims,dims,filespace,hdf_error)

      data_count(1)=n2
      data_count(2)=kend-kstart+1

      data_offset(1)=0
      data_offset(2)=kstart-1

     !pressure
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id,hdf_error)

      CALL h5pset_fapl_mpio_f(plist_id, comm, info,hdf_error)

      CALL h5fcreate_f(filnam5, H5F_ACC_TRUNC_F, file_id, &
         hdf_error, access_prp=plist_id)

      CALL h5pclose_f(plist_id, hdf_error)

      CALL h5dcreate_f(file_id, 'pr', H5T_NATIVE_DOUBLE, &
          filespace,dset_pr, hdf_error)

      CALL h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      CALL h5dget_space_f(dset_pr, slabspace, hdf_error)
      
      CALL h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                  data_offset, data_count, hdf_error)
      
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      
      CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)
      
      CALL h5dwrite_f(dset_pr, H5T_NATIVE_DOUBLE, &
         pr(1:n2,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, &
         mem_space_id = memspace, xfer_prp = plist_id)
      
      CALL h5pclose_f(plist_id, hdf_error)

      CALL h5dclose_f(dset_pr, hdf_error)

      CALL h5sclose_f(memspace, hdf_error)
      
      CALL h5fclose_f(file_id, hdf_error)

    
     !q2
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)

      CALL h5pset_fapl_mpio_f(plist_id, comm, info,  hdf_error)

      CALL h5fcreate_f(filnam3, H5F_ACC_TRUNC_F, file_id, &
           hdf_error, access_prp=plist_id)

      CALL h5pclose_f(plist_id, hdf_error)

      CALL h5dcreate_f(file_id, 'Vy', H5T_NATIVE_DOUBLE, &
                     filespace, dset_q2, hdf_error)

      CALL h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      CALL h5dget_space_f(dset_q2, slabspace, hdf_error)

      CALL h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                        data_offset, data_count, hdf_error)

      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                         hdf_error)
      CALL h5dwrite_f(dset_q2, H5T_NATIVE_DOUBLE, &
        q2(1:n2,kstart:kend), dims,  &
        hdf_error, file_space_id = slabspace, &
        mem_space_id = memspace, xfer_prp = plist_id)

      CALL h5pclose_f(plist_id, hdf_error)

      CALL h5dclose_f(dset_q2, hdf_error)

      CALL h5sclose_f(memspace, hdf_error)
      
      CALL h5fclose_f(file_id, hdf_error)


     !q3
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)

      CALL h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

      CALL h5fcreate_f(filnam4, H5F_ACC_TRUNC_F, file_id, &
           hdf_error, access_prp=plist_id)

      CALL h5pclose_f(plist_id, hdf_error)

      CALL h5dcreate_f(file_id, 'Vz', H5T_NATIVE_DOUBLE, &
                  filespace, dset_q3, hdf_error)

      CALL h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      CALL h5dget_space_f(dset_q3, slabspace, hdf_error)
      
      CALL h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                     data_offset, data_count, hdf_error)
      
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                     hdf_error)
      
      CALL h5dwrite_f(dset_q3, H5T_NATIVE_DOUBLE, &
        q3(1:n2,kstart:kend), dims,  &
        hdf_error, file_space_id = slabspace, &
        mem_space_id = memspace, xfer_prp = plist_id)

      CALL h5pclose_f(plist_id, hdf_error)

      CALL h5dclose_f(dset_q3, hdf_error)

      CALL h5sclose_f(memspace, hdf_error)
      
      CALL h5fclose_f(file_id, hdf_error)
      
     
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

     
     !ru2
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)

      CALL h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

      CALL h5fcreate_f(filnam7, H5F_ACC_TRUNC_F, file_id, &
           hdf_error, access_prp=plist_id)

      CALL h5pclose_f(plist_id, hdf_error)

      CALL h5dcreate_f(file_id, 'ru2', H5T_NATIVE_DOUBLE, &
                  filespace, dset_q3, hdf_error)

      CALL h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      CALL h5dget_space_f(dset_q3, slabspace, hdf_error)
      
      CALL h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                     data_offset, data_count, hdf_error)
      
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                     hdf_error)
      
      CALL h5dwrite_f(dset_q3, H5T_NATIVE_DOUBLE, &
        ru2(1:n2,kstart:kend), dims,  &
        hdf_error, file_space_id = slabspace, &
        mem_space_id = memspace, xfer_prp = plist_id)

      CALL h5pclose_f(plist_id, hdf_error)

      CALL h5dclose_f(dset_q3, hdf_error)

      CALL h5sclose_f(memspace, hdf_error)
      
      CALL h5fclose_f(file_id, hdf_error)
     
     
     !ru3
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)

      CALL h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

      CALL h5fcreate_f(filnam8, H5F_ACC_TRUNC_F, file_id, &
           hdf_error, access_prp=plist_id)

      CALL h5pclose_f(plist_id, hdf_error)

      CALL h5dcreate_f(file_id, 'ru3', H5T_NATIVE_DOUBLE, &
                  filespace, dset_q3, hdf_error)

      CALL h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      CALL h5dget_space_f(dset_q3, slabspace, hdf_error)
      
      CALL h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                     data_offset, data_count, hdf_error)
      
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                     hdf_error)
      
      CALL h5dwrite_f(dset_q3, H5T_NATIVE_DOUBLE, &
        ru3(1:n2,kstart:kend), dims,  &
        hdf_error, file_space_id = slabspace, &
        mem_space_id = memspace, xfer_prp = plist_id)

      CALL h5pclose_f(plist_id, hdf_error)

      CALL h5dclose_f(dset_q3, hdf_error)

      CALL h5sclose_f(memspace, hdf_error)
      
      CALL h5fclose_f(file_id, hdf_error)
       

     !write some scalars for restart consistency and grid nodes
      if (myid .eq. 0) then
        
       
       !write grid files
        filnamgrid = trim(outfol)//'/restart/ycoord.dat'
        
        open(151,file=filnamgrid,status='unknown')
        do i=1, n2
          write(151,*) rc(i)
        enddo
        close(151)
        
        filnamgrid =trim(outfol)//'/restart/zcoord.dat'
        
        open(151,file=filnamgrid,status='unknown')
        do i=1, n3
          write(151,*) zz(i)
        enddo
        close(151)
       
        print*,''
        print*, ' ... written restart files'
        print*,''
  
      endif
      
      
      ENDSUBROUTINE
      
     !===================================================================

     !MPI read arrays with ghost cells from .h5 
      
      SUBROUTINE mpi_read_continua(n2o,n3o,ks,ke,intvar,qua)
      
      USE mpih
      USE param
      USE hdf5
      
      implicit none
     !---------------------------------------------------------- 
      integer, intent(in) :: ks,ke,n2o,n3o
      real, dimension(1:n2o,ks-lvlhalo:ke+lvlhalo)::qua

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_qua

      integer(HSIZE_T) :: dims(2)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(2) :: data_count  
      integer(HSSIZE_T), dimension(2) :: data_offset 

      integer :: comm, info
      integer :: ndims

      integer, intent(in) :: intvar
      character*70 :: filnam1
      character*10 :: dsetname
     !----------------------------------------------------------


      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

     !Select file and dataset based on intvar
      selectcase (intvar)
        case(2)
          dsetname = trim('Vy')
          filnam1 = 'input_FSI/restart/rest_q2.h5'
        case(3)
          dsetname = trim('Vz')
          filnam1 = 'input_FSI/restart/rest_q3.h5'
        case(4)
          dsetname = trim('pr')
          filnam1 = 'input_FSI/restart/rest_pr.h5'
        case(5)  
          dsetname = trim('dens')
          filnam1 = 'input_FSI/restart/rest_dens.h5'
      endselect

     !Set offsets and element counts
      ndims=2

      dims(1)=n2o
      dims(2)=n3o-1

      data_count(1)=n2o
      data_count(2)=ke-ks+1

      data_offset(1)=0
      data_offset(2)=ks-1

        
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, &
            hdf_error)


      CALL h5pset_fapl_mpio_f(plist_id, comm, info, &
            hdf_error)

      CALL h5fopen_f(filnam1, H5F_ACC_RDONLY_F, file_id, &
            hdf_error, access_prp=plist_id)

      CALL h5dopen_f(file_id, dsetname,dset_qua, hdf_error)

      CALL h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      CALL h5dget_space_f(dset_qua, slabspace, hdf_error)
      
      CALL h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                      data_offset, data_count, hdf_error)
      
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      
      CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)
      
      CALL h5dread_f(dset_qua, H5T_NATIVE_DOUBLE, &
        qua(1:n2o,ks:ke), dims,  &
        hdf_error, file_space_id = slabspace, &
        mem_space_id = memspace, xfer_prp = plist_id)
      
      CALL h5pclose_f(plist_id, hdf_error)

      CALL h5dclose_f(dset_qua, hdf_error)

      CALL h5sclose_f(memspace, hdf_error)
      
      CALL h5fclose_f(file_id, hdf_error)


      ENDSUBROUTINE

     !===================================================================

     !MPI read arrays withoout ghost cells from .h5 

      SUBROUTINE mpi_read_continua2(n2o,n3o,ks,ke,intvar,qua)
      
      USE mpih
      USE param
      USE hdf5
      
      implicit none
     !---------------------------------------------------------- 
      integer, intent(in) :: ks,ke,n2o,n3o
      real, dimension(1:n2o,ks:ke)::qua

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_qua

      integer(HSIZE_T) :: dims(2)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(2) :: data_count  
      integer(HSSIZE_T), dimension(2) :: data_offset 

      integer :: comm, info
      integer :: ndims

      integer, intent(in) :: intvar
      character*70 :: filnam1
      character*10 :: dsetname
     !----------------------------------------------------------


      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

     !Select file and dataset based on intvar
      selectcase (intvar)
        case(6)  
          dsetname = trim('ru1')
          filnam1 = 'input_FSI/restart/rest_ru1.h5'
        case(7)  
          dsetname = trim('ru2')
          filnam1 = 'input_FSI/restart/rest_ru2.h5'
        case(8)  
          dsetname = trim('ru3')
          filnam1 = 'input_FSI/restart/rest_ru3.h5'
      endselect

     !Set offsets and element counts
      ndims=2

      dims(1)=n2o
      dims(2)=n3o-1

      data_count(1)=n2o
      data_count(2)=ke-ks+1

      data_offset(1)=0
      data_offset(2)=ks-1

        
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, &
            hdf_error)


      CALL h5pset_fapl_mpio_f(plist_id, comm, info, &
            hdf_error)

      CALL h5fopen_f(filnam1, H5F_ACC_RDONLY_F, file_id, &
            hdf_error, access_prp=plist_id)

      CALL h5dopen_f(file_id, dsetname,dset_qua, hdf_error)

      CALL h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      CALL h5dget_space_f(dset_qua, slabspace, hdf_error)
      
      CALL h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                      data_offset, data_count, hdf_error)
      
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      
      CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)
      
      CALL h5dread_f(dset_qua, H5T_NATIVE_DOUBLE, &
        qua(1:n2o,ks:ke), dims,  &
        hdf_error, file_space_id = slabspace, &
        mem_space_id = memspace, xfer_prp = plist_id)
      
      CALL h5pclose_f(plist_id, hdf_error)

      CALL h5dclose_f(dset_qua, hdf_error)

      CALL h5sclose_f(memspace, hdf_error)
      
      CALL h5fclose_f(file_id, hdf_error)

      ENDSUBROUTINE

     !===================================================================




