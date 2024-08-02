      
     !It write on -h5 files the nodes grid and the cell-center grid

      SUBROUTINE inimov
      
      USE mpih
      USE param
      USE hdf5

      implicit none
     !----------------------------------------------------------
      integer hdf_error
      integer(HID_T) :: file_grid
      integer(HID_T) :: dset_grid
      integer(HID_T) :: dspace_grid
      integer(HSIZE_T) :: dims_grid(1)
      character*70 :: namfile_c, namfile_n, namfol
      character(50) :: filen
      logical :: dir_e
     !----------------------------------------------------------


      if (myid.eq.0) then

        namfile_c=trim(outfol)//'/movie/cell_center_grid.h5'
        namfile_n=trim(outfol)//'/movie/nodes_grid.h5' 
  
        namfol=trim(outfol)//'/movie'
        
        inquire(file=trim(namfol), exist=dir_e)
        if (.not.dir_e) then
          CALL system('mkdir '//trim(namfol))
        endif


        filen=trim(outfol)//'/movie/movierb.out'
        open(319, file=trim(filen), action="write",status='unknown')
  

       !write the cell center grid
  
        CALL h5fcreate_f(namfile_c,H5F_ACC_TRUNC_F,file_grid,hdf_error)
  
        dims_grid(1)=n2m
        CALL h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
  
        CALL h5dcreate_f(file_grid, 'rm', H5T_NATIVE_DOUBLE, &
                  dspace_grid, dset_grid, hdf_error)
  
        CALL h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, rm(1:n2m),    & 
                  dims_grid,hdf_error)
  
        CALL h5dclose_f(dset_grid, hdf_error)
        CALL h5sclose_f(dspace_grid, hdf_error)
  
        dims_grid(1)=n3m
        CALL h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
        CALL h5dcreate_f(file_grid, 'zm', H5T_NATIVE_DOUBLE,   &    
                  dspace_grid, dset_grid, hdf_error)
  
        CALL h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, zm(1:n3m),   &
                  dims_grid, hdf_error)
  
  
        CALL h5dclose_f(dset_grid, hdf_error)
        CALL h5sclose_f(dspace_grid, hdf_error)
  
        CALL h5fclose_f(file_grid, hdf_error)
        
  
       !write the nodes grid
  
        CALL h5fcreate_f(namfile_n,H5F_ACC_TRUNC_F,file_grid,hdf_error)
  
        dims_grid(1)=n2
        CALL h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
  
        CALL h5dcreate_f(file_grid, 'rc', H5T_NATIVE_DOUBLE,  & 
                  dspace_grid, dset_grid, hdf_error)
  
        CALL h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, rc(1:n2),  &
                  dims_grid,hdf_error)
  
        CALL h5dclose_f(dset_grid, hdf_error)
        CALL h5sclose_f(dspace_grid, hdf_error)
  
        dims_grid(1)=n3
        CALL h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
        CALL h5dcreate_f(file_grid, 'zz', H5T_NATIVE_DOUBLE,     &
                  dspace_grid, dset_grid, hdf_error)
  
        CALL h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, zz(1:n3),    &
                  dims_grid, hdf_error)
  
  
        CALL h5dclose_f(dset_grid, hdf_error)
        CALL h5sclose_f(dspace_grid, hdf_error)
  
        CALL h5fclose_f(file_grid, hdf_error)

      endif

      ENDSUBROUTINE


