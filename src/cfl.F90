      
     !It computes the maximum local CFL number 

      SUBROUTINE cfl

      USE param
      USE local_arrays,  ONLY: q2,q3
      USE mpih
      USE mpi_param,  ONLY: kstart,kend
      
      implicit none
     !----------------------------------------------------------
      real    :: my_cflm, myq3avg
      integer :: j, k, jp, kp, i, ip
      real :: qcf, udx3, udx2
     !----------------------------------------------------------


      my_cflm=1.0e-8 
      myq3avg=0.0

                                                                       
      do k=kstart, kend
        udx3=udx3m(k)
        kp=k+1

        do j=1, n2m
          udx2=udx2m(j)
          jp=j+1
          
          qcf=( abs((q2(j,k)+q2(jp,k))*0.5*udx2)+  &
                abs((q3(j,k)+q3(j,kp))*0.5*udx3))
          my_cflm = max(my_cflm,qcf)
          myq3avg = myq3avg+q3(j,k)*(zz(kp)-zz(k))*  &
                    (rc(jp)-rc(j))

        enddo
      enddo

            
      CALL MPI_ALLREDUCE(my_cflm,cflm,1,MDP,MPI_MAX,  &
               MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(myq3avg,q3avg,1,MDP,MPI_SUM, &
               MPI_COMM_WORLD,ierr)


     !volume averaged value
      q3avg=q3avg/(rext2*alx3)

       
      ENDSUBROUTINE                        

