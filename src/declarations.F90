
      MODULE param
        
        implicit none
        
       !input variables (from fluid_par.in, global_par.in)
        integer   :: n2, n3
        integer   :: nsst
        real      :: tprint
        real      :: alx3
        real      :: rext, rext2
        real      :: dt, dt_o, resid, cflmax, tsta
        real      :: dtmax, cfllim 
        integer   :: nson
        integer   :: dimSD
        real      :: sclf, rhop
        real      :: q3uw, q3lw 
        integer   :: buoy
        real      :: cou

       !new variables and redefined variables 
        logical :: idtv, nwrit, stopsim
        logical :: fcut
        logical :: wcheck
        logical :: starea
        logical :: inslws,inslwn
        integer :: wfun
        real :: cflm, mflr
        real :: fortime
        integer :: shellG
       
       !grid parameters
        integer :: n2m, n3m
        integer :: n2mh, n3mh
        integer :: n2mp, n3mp
        real :: dx2,dx3 !,dx1
        real :: dx2q,dx3q,dx1q
        real :: str2, str3, n3strm,n3strp, n2strm,n2strp
        integer :: istr3, istr2, n3str, n2str
        real :: str3u, str3l, str2u, str2l
        real, dimension(2) :: int3, int2
        
!        real, dimension(:), allocatable :: tc,tm
        real, dimension(:), allocatable :: rc,rm,g2rc,g2rm,etay,etaym
        real, dimension(:), allocatable :: zz,zm,g3rc,g3rm,etaz,etazm
        
       !quantities for derivatives 
        real, dimension(:), allocatable :: udx3c, udx3m 
        real, dimension(:), allocatable :: udx2c, udx2m

       !grid indices
        integer, dimension(:), allocatable :: jmv,jpv,jpc,jmc
!        integer, dimension(:), allocatable :: imv,ipv
        integer, dimension(:), allocatable :: imhv, jmhv
        integer, dimension(:), allocatable :: kmc,kpc,kmv,kpv,kup,kum
       
       !metric coefficients 
        real, dimension(:), allocatable :: ap2j,ac2j,am2j
        real, dimension(:), allocatable :: ap3j,ac3j,am3j
        real, dimension(:), allocatable :: apscj,acscj,amscj
        real, dimension(:), allocatable :: ap3ck,ac3ck,am3ck
        real, dimension(:), allocatable :: ap3sk,ac3sk,am3sk
        real, dimension(:), allocatable :: ap3ssk,ac3ssk,am3ssk

       !variables for FFTW and Poisson solver
        integer :: fwd_plan,bck_plan
        real, dimension(13) :: ifx1
        
        real, dimension(:), allocatable :: trigx1
        real, dimension(:), allocatable :: ak2,ap
        real, dimension(:), allocatable :: ak1,ao
        real, dimension(:), allocatable :: amphk,acphk,apphk
        real, dimension(:), allocatable :: amphj,apphj,acphj
       
       !other variables 
        integer  :: iaxsy
        real :: rint, q3avg
        real :: pi, ggg
        real :: al,ga,ro 
        real :: beta
        real :: qqmax,qqtot
        real :: re
        real :: anusslow,anussupp
        real :: denmax,denmin,densm
        real, dimension(1:3) :: gam,rom,alm
        real, dimension(:), allocatable :: denbs,denbn
        integer, dimension(:), allocatable :: iwnt,iwst
        real, dimension(50) :: nors
        
       !nondimensional number and scale parameter
        real :: ren, visc

       !time variables
        real :: time
        real :: ntime, ntime1, ntime2
        real :: trestart
        real, dimension(2) :: tit_f    !time for one fluid time step
        real, dimension(2) :: tit_s    !time for one shell time step
        real, dimension(2) :: ttot     !time for simulation
        real, dimension(2) :: tload, tins

       !output folders and output parameters
        character(50) :: dirout_f,dirout_s
        character(50) :: outfol
        real :: tframes, tframef, tcond 

       !relevant flow control variables
        logical :: ibact, irest

       !for Lapack
        real, dimension(:), allocatable :: WR 
        real, dimension(:,:), allocatable :: zmm,zmmt
        
              
      ENDMODULE param

     !====================================================================                

      MODULE outflow_vars

      USE param

        implicit none

        real :: qinf, qout
!        real, dimension(:,:), allocatable :: qb1s, qb1n, dqb1s, dqb1n
        real, dimension(:), allocatable :: qb2s, qb2n, dqb2s, dqb2n
        real, dimension(:), allocatable :: qb3s, qb3n, dqb3s, dqb3n
        real, dimension(:), allocatable :: dq2x2o, dq3x2o
        real, dimension(:), allocatable :: dq2x2on, dq3x2on
        real, dimension(:), allocatable :: dq2x2os, dq3x2os

       !for strong coupling
        real, dimension(:), allocatable :: qb2s_o, qb3s_o
        real, dimension(:), allocatable :: qb2n_o, qb3n_o
        real, dimension(:), allocatable :: dq2x2o_o, dq3x2o_o

      ENDMODULE

     !====================================================================

      MODULE local_arrays
      
      USE param
        
        implicit none
        
        real, allocatable, dimension(:,:) :: q2, q3
        real, allocatable, dimension(:,:) :: pr, dens, hro, rhs
        real, allocatable, dimension(:,:) :: ru2, ru3, ruro
        real, allocatable, dimension(:,:) :: qcap
        real, allocatable, dimension(:,:) :: dph
        real, allocatable, dimension(:,:) :: rhs2, rhs3
        real, allocatable, dimension(:,:) :: q2f, q3f

       !for strong coupling
        real, allocatable, dimension(:,:) :: q1_o, q2_o, q3_o
        real, allocatable, dimension(:,:) :: pr_o, dens_o
        real, allocatable, dimension(:,:) :: ru1_o, ru2_o, ru3_o
      
      ENDMODULE local_arrays
     
!     !====================================================================                
      
      MODULE mpih

        USE mpi
       
        implicit none
        
!        INCLUDE 'mpif.h'
        
        integer :: myid, numtasks, ierr
        integer, parameter :: master=0

       !three ghost cell rows are needed for the support domain of 
       !MLS forcing
        integer :: lvlhalo  
        integer :: MDP = MPI_DOUBLE_PRECISION
        integer :: STATUS(MPI_STATUS_SIZE,4)
        integer :: req(1:4)
        integer(kind=MPI_OFFSET_KIND) :: disp, offset
      
      ENDMODULE mpih

     !==================================================================== 

      MODULE mpi_param
        
        implicit none
        
        integer :: istart, iend, jstart, jend, kstart, kend
        integer :: jstartp,jendp
        integer :: di, dj, dk, mydata, mydatam
        integer :: djp
        integer, allocatable, dimension(:) :: offseti, offsetj, offsetk
        integer, allocatable, dimension(:) :: offsetjp
        integer, allocatable, dimension(:) :: counti, countj, countk
        integer, allocatable, dimension(:) :: countjp
        integer, allocatable, dimension(:) :: countf
        integer, allocatable, dimension(:) :: offsetf 
      
      ENDMODULE mpi_param

     !==================================================================== 
     
      MODULE mls_param
        
      implicit none

      integer :: lgrdis
      integer :: forst, nlgrtot
      integer :: bfo
      real :: radiusSD, kprobe, ddmin
      real :: tprtag
      real :: scaledmin
      real :: nbcu, nbcv, nbcw
      
      logical :: cbody,prlgr

      integer, allocatable :: eulins(:,:,:,:), eulin_L(:,:,:)
      real, allocatable :: cl(:,:), phi(:,:,:), phipred(:,:,:)

      real, allocatable :: phi_L(:,:,:)
      real, allocatable :: press(:), press_r(:)
      real, allocatable :: tau(:,:), tau_r(:,:)

      real, dimension(:,:), allocatable :: fey, fez
      real, dimension(:), allocatable :: bcerru, bcerrv, bcerrw
      real, dimension(:), allocatable :: bcerru_r, bcerrv_r, bcerrw_r
      real, dimension(:), allocatable :: bcerru_tot,bcerrv_tot,bcerrw_tot

      integer, dimension(:), allocatable :: ftag1, ftag2, ftag3
      integer, dimension(:), allocatable :: ftag_g1, ftag_g2, ftag_g3
        
      ENDMODULE mls_param

     !==================================================================== 
      
      MODULE rbm

      real :: rbden, rba(2), ideti, rrr
      real :: rbgravF(2), forgcm(2)

      integer :: rbmcon(3)

      real, dimension(3) :: ipert
      real, dimension(3) :: cmp_np1, cmv_np1, cma_np1  
      real, dimension(3) :: cmp_n, cmv_n, cma_n 
      real, dimension(3) :: cmp_k, cmv_k, cma_k, cma_km1, cma_p, cma_c
      real, dimension(3) :: cma_nm1
      real, dimension(3) :: cmp_np1_norel
      real, dimension(3) :: kk,cc,kk3
      real, dimension(2) :: cmposi 
      real :: cmroti,cmrot,cmtra(2)
      real :: momgcm,momcm
      real :: kinori(3,3)

      integer :: nlm
      integer, allocatable :: lmind(:,:,:)
      real, allocatable :: lmcc(:,:),lmA(:),lmnor(:,:),vert(:,:)
      real, allocatable :: lmvel(:,:),lmacc(:,:)
      real, dimension(:,:), allocatable :: Ftot

      logical :: rtcoup 

      ENDMODULE
      
     !==================================================================== 

      MODULE fsicoup

      integer :: fsic       
      real :: atk, sceps
      real, dimension(:), allocatable :: dsc_k, vsc_k, asc_k, fsc_k
      real, dimension(:), allocatable :: asc_c, asc_km1, asc_p

       ENDMODULE

     !==================================================================== 


