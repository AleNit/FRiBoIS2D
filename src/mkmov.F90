           
     !It writes to file primitive variables on a domain slice 
     !with tframef intervals
      
      SUBROUTINE mkmov_slice(ns)

      USE local_arrays  
      USE mpi_param
      USE mpih
      USE hdf5
      USE param

      IMPLICIT NONE
     !----------------------------------------------------------      
      integer, intent(in) :: ns 
     !---------------------------------------------------------- 
      integer :: ic, jc, kc, ip, jp, kp
      integer :: i1, i2, j1, j2, k1, k2
      integer :: ndims, itime, comm, info
      integer :: hdf_error
      real :: tprfi
      integer(HID_T) :: filespace, memspace, plist_id, file_id
      integer(HID_T) :: dset_q1v, dset_q2v, dset_q3v, dset_prv
      integer(HSIZE_T) :: dims(2), dimsp(2)
      integer(HSIZE_T), dimension(2) :: data_count, data_countp
      integer(HSSIZE_T), dimension(2) :: data_offset
      character(70) :: filnam1,filnamxdm
      character(5) :: ipfi
     !----------------------------------------------------------      


     !Sort out MPI definitions and file names
      tprfi = 1/tframef
      itime=nint(time*tprfi)
      write(ipfi,82)itime
   82 format(i5.5)


      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

      ndims=2


      filnam1=trim(outfol)//'/movie/frame_'//ipfi//'.h5'

      dims=[n2,n3]
      dimsp=[n2m,n3m]
      data_count=[n2,kend-kstart+1]
      data_countp=[n2m,kend-kstart+1]
      data_offset=[0,kstart-1]


     !Open file and create dataspace
      CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdf_error)
      CALL h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      CALL h5fcreate_f(filnam1, H5F_ACC_TRUNC_F, file_id, hdf_error, & 
                 access_prp=plist_id)
      CALL h5pclose_f(plist_id, hdf_error)
      
      CALL h5screate_simple_f(ndims, dims,filespace, hdf_error)


     !------------------------------------------- velocity components
      CALL h5dcreate_f(file_id,'Vy',H5T_NATIVE_DOUBLE,filespace,    &
                      dset_q2v, hdf_error) 
      CALL h5dcreate_f(file_id,'Vz',H5T_NATIVE_DOUBLE,filespace,    &
                      dset_q3v, hdf_error) 

      CALL h5screate_simple_f(ndims, data_count, memspace, hdf_error)


      CALL h5dget_space_f(dset_q2v, filespace, hdf_error)
      CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,     &
                            data_offset, data_count, hdf_error)
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,  &
                              hdf_error)

      
      CALL h5dwrite_f(dset_q2v, H5T_NATIVE_DOUBLE,   &
          q2(1:n2,kstart:kend), dims,    &
          hdf_error, file_space_id = filespace, mem_space_id = memspace, &
          xfer_prp = plist_id)
      CALL h5pclose_f(plist_id, hdf_error)

 

      CALL h5dget_space_f(dset_q3v, filespace, hdf_error)
      CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,   &
                            data_offset, data_count, hdf_error)
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,  &
                              hdf_error)
      CALL h5dwrite_f(dset_q3v, H5T_NATIVE_DOUBLE,  &
          q3(1:n2,kstart:kend), dims, &
          hdf_error, file_space_id = filespace, mem_space_id = memspace, &
          xfer_prp = plist_id)
      CALL h5pclose_f(plist_id, hdf_error)
 
      CALL h5dclose_f(dset_q2v, hdf_error)
      CALL h5dclose_f(dset_q3v, hdf_error)



     !------------------------------------------- pressure
      CALL h5screate_simple_f(ndims, dimsp,filespace, hdf_error)


     !Create the dataset with default properties
      CALL h5dcreate_f(file_id,'Pr',H5T_NATIVE_DOUBLE,filespace,    &
                      dset_prv, hdf_error)

      CALL h5screate_simple_f(ndims, data_countp, memspace, hdf_error)
     

      CALL h5dget_space_f(dset_prv, filespace, hdf_error)
      CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,   &
                            data_offset, data_countp, hdf_error)
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,  &
                             hdf_error)
      
      CALL h5dwrite_f(dset_prv, H5T_NATIVE_DOUBLE,  &
          pr(1:n2m,kstart:kend), dimsp,   &
          hdf_error, file_space_id = filespace, mem_space_id = memspace, &
          xfer_prp = plist_id)
      CALL h5pclose_f(plist_id, hdf_error)

 
      CALL h5dclose_f(dset_prv, hdf_error)
 
      CALL h5sclose_f(memspace, hdf_error)
      CALL h5sclose_f(filespace, hdf_error)
      CALL h5fclose_f(file_id, hdf_error)


      ENDSUBROUTINE                     

