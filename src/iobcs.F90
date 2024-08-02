
     !It sets inflow/outflow boundary conditions through 
     !wave equation

      SUBROUTINE boucdq
      
      USE param
      USE local_arrays,  only: q2, q3, pr
      USE mpi_param,    only: kstart, kend
      USE mpih
      USE outflow_vars
      
      implicit none
     !---------------------------------------------------------- 
      integer :: i,j
      real :: dq1x2, dq2x2, dq3x2
      real :: area, darea
      real :: qout2, qoutf, cor, velin, delta
     !---------------------------------------------------------- 


      dqb2s(:)=0.0


     !steady inflow condition 
      dqb3s(:)=0.0     


     !Radiative B.C. at outflow 
     !The spatial derivative is updated in time following the 
     !convective term time scheme
      if (kend.eq.n3m) then
        do j=1, n2m
          
          dq2x2=(qb2n(j)-q2(j,n3m))*2.0*udx3c(n3)
          dqb2n(j)=-dt*(ga*dq2x2+ro*dq2x2o(j))*cou

          dq3x2=(qb3n(j)-q3(j,n3m))*udx3c(n3)
          dqb3n(j)=-dt*(ga*dq3x2+ro*dq3x2o(j))*cou

         !store spatial derivative for next time step
          dq2x2o(j)=dq2x2
          dq3x2o(j)=dq3x2
       
        enddo   
      endif

      CALL MPI_BCAST(dqb2n,n2,MDP,numtasks-1,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dqb3n,n2,MDP,numtasks-1,MPI_COMM_WORLD,ierr)


     !Flow is adjusted so that it is globally divergence free
     !First mass imbalance is calculated
      area=0.0
      qout=0.0
      qinf=0.0
      do j=1,n2m
        darea=g2rm(j)/dx2
        area = area+darea
        qout=qout+dqb3n(j)*darea
        qinf=qinf+dqb3s(j)*darea
      enddo

     !Then the correction is calculated and implemented at the outflow
      cor=(qinf-qout)/area
      qout2=0.0
      qoutf=0.0
      do j=1,n2m
        darea=g2rm(j)/dx2
        area = area+darea
        dqb3n(j)=dqb3n(j)+cor
        qout2=qout2+cor*darea
        qoutf=qoutf+dqb3n(j)*darea
      enddo

      
      ENDSUBROUTINE

     !===================================================================
      
      SUBROUTINE boucqt
      
      USE param
      USE outflow_vars
      USE mpih
      
      implicit none
     !---------------------------------------------------------- 
      real :: darea
      integer :: i, j      
     !---------------------------------------------------------- 

      qout=0.0
      qinf=0.0

     !Inflow/outflow velocities are updated, and mass imbalance is checked
      do j=1,n2m
        darea=g2rm(j)/dx2
        qb2s(j) = qb2s(j)+dqb2s(j)
        qb3s(j) = qb3s(j)+dqb3s(j)
        qb2n(j) = qb2n(j)+dqb2n(j)
        qb3n(j) = qb3n(j)+dqb3n(j)
        qout=qout+qb3n(j)*darea
        qinf=qinf+qb3s(j)*darea
      enddo

      
      ENDSUBROUTINE
     
     !===================================================================
 
     !write inflow/outflow arrays to file for the restart procedure   

      SUBROUTINE writeqbns
      
      USE param
      USE outflow_vars
      USE mpi_param,    only: kstart, kend
      USE mpih
      USE hdf5
      
      implicit none
     !---------------------------------------------------------- 
      integer :: hdf_error, ndims
      character(100) :: filename
      integer(HSIZE_T) :: dims(1)
      integer(HID_T) :: filespace, memspace, dset
     !---------------------------------------------------------- 

      
      if (kend==n3m) then

        filename=trim(outfol)//'/restart/outflow.h5'
  
        ndims=1
        dims=[n2]
        
  
        CALL h5fcreate_f(filename,H5F_ACC_TRUNC_F,filespace,hdf_error)
  
        CALL h5screate_simple_f(ndims,dims,memspace, hdf_error)
  
 
        CALL h5dcreate_f(filespace,'qb2n',H5T_NATIVE_DOUBLE,  &
                         memspace,dset,hdf_error)
        CALL h5dwrite_f(dset,H5T_NATIVE_DOUBLE,qb2n,dims,hdf_error)
        CALL h5dclose_f(dset, hdf_error)
       
        CALL h5dcreate_f(filespace,'qb3n',H5T_NATIVE_DOUBLE,  &
                         memspace,dset,hdf_error)
        CALL h5dwrite_f(dset,H5T_NATIVE_DOUBLE,qb3n,dims,hdf_error)
        CALL h5dclose_f(dset, hdf_error)
  
 
        CALL h5dcreate_f(filespace,'dq2x2o',H5T_NATIVE_DOUBLE,  &
                         memspace,dset,hdf_error)
        CALL h5dwrite_f(dset,H5T_NATIVE_DOUBLE,dq2x2o,dims,hdf_error)
        CALL h5dclose_f(dset, hdf_error)
        
        CALL h5dcreate_f(filespace,'dq3x2o',H5T_NATIVE_DOUBLE,  &
                         memspace,dset,hdf_error)
        CALL h5dwrite_f(dset,H5T_NATIVE_DOUBLE,dq3x2o,dims,hdf_error)
        CALL h5dclose_f(dset, hdf_error)
  
  
        CALL h5sclose_f(memspace, hdf_error)
        
        CALL h5fclose_f(filespace, hdf_error)

      endif


      if (kstart==1) then
      
        filename=trim(outfol)//'/restart/inflow.h5'
  
        ndims=1
        dims=[n2]
        
  
        CALL h5fcreate_f(filename,H5F_ACC_TRUNC_F,filespace,hdf_error)
  
        CALL h5screate_simple_f(ndims,dims,memspace, hdf_error)
  
        CALL h5dcreate_f(filespace,'qb2s',H5T_NATIVE_DOUBLE,  &
                         memspace,dset,hdf_error)
        CALL h5dwrite_f(dset,H5T_NATIVE_DOUBLE,qb2s,dims,hdf_error)
        CALL h5dclose_f(dset, hdf_error)
       
        CALL h5dcreate_f(filespace,'qb3s',H5T_NATIVE_DOUBLE,  &
                         memspace,dset,hdf_error)
        CALL h5dwrite_f(dset,H5T_NATIVE_DOUBLE,qb3s,dims,hdf_error)
        CALL h5dclose_f(dset, hdf_error)
        
        CALL h5sclose_f(memspace, hdf_error)
        
        CALL h5fclose_f(filespace, hdf_error)

      endif
      
      
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

     
      ENDSUBROUTINE
     
     !===================================================================
      
      SUBROUTINE readqbns
      
      USE param
      USE outflow_vars
      USE mpih
      USE mpi_param,  ONLY: kstart, kend
      USE hdf5
      
      implicit none
     !---------------------------------------------------------- 
      integer :: hdf_error, ndims
      character(100) :: filename
      integer(HSIZE_T) :: dims(1)
      integer(HID_T) :: dspace_id, dset_id, file_id
     !---------------------------------------------------------- 

      
      if (kstart==1) then
      
        filename='input_FSI/restart/inflow.h5'
  
        ndims=1
        dims=[n2]
  
        CALL h5open_f(hdf_error) 
      
        CALL h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf_error)
    
        CALL h5dopen_f(file_id,'qb2s',dset_id,hdf_error)
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,qb2s,dims,hdf_error)
        CALL h5dclose_f(dset_id,hdf_error)
  
        CALL h5dopen_f(file_id,'qb3s',dset_id,hdf_error)
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,qb3s,dims,hdf_error)
        CALL h5dclose_f(dset_id,hdf_error)
  
        CALL h5fclose_f(file_id,hdf_error)

      endif


      if (kend==n3m) then

        filename='input_FSI/restart/outflow.h5'
  
        ndims=1
        dims=[n2]
  
        CALL h5open_f(hdf_error) 
      
        CALL h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf_error)
      
        CALL h5dopen_f(file_id,'qb2n',dset_id,hdf_error)
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,qb2n,dims,hdf_error)
        CALL h5dclose_f(dset_id,hdf_error)
  
        CALL h5dopen_f(file_id,'qb3n',dset_id,hdf_error)
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,qb3n,dims,hdf_error)
        CALL h5dclose_f(dset_id,hdf_error)
        
        CALL h5dopen_f(file_id,'dq2x2o',dset_id,hdf_error)
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,dq2x2o,dims,hdf_error)
        CALL h5dclose_f(dset_id,hdf_error)
  
        CALL h5dopen_f(file_id,'dq3x2o',dset_id,hdf_error)
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,dq3x2o,dims,hdf_error)
        CALL h5dclose_f(dset_id,hdf_error)
  
        CALL h5fclose_f(file_id,hdf_error)
 
      endif


      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)


      CALL MPI_BCAST(qb2s,n2,MDP,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(qb3s,n2,MDP,0,MPI_COMM_WORLD,ierr)
      
      CALL MPI_BCAST(qb2n,n2,MDP,numtasks-1,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(qb3n,n2,MDP,numtasks-1,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dq2x2o,n2,MDP,numtasks-1,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dq3x2o,n2,MDP,numtasks-1,MPI_COMM_WORLD,ierr)

     
      ENDSUBROUTINE

