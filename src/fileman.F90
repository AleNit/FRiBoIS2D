      
     !It open output files and create necessary output folders
      
      SUBROUTINE openfi

      USE param
      USE mpih,  ONLY: numtasks, myid
      USE mls_param,  ONLY: cbody,prlgr
      
      implicit none
     !---------------------------------------------------------- 
      integer :: i
      logical :: dir_e
      integer, dimension(8) :: det
      character(10), dimension(3) :: b
      character(50) :: filen, newdir
      character(1) :: patch_ID
     !---------------------------------------------------------- 

     
      filen=trim(outfol)//'/restart'
      inquire(file=filen, exist=dir_e)
      if (.not.dir_e) then
        CALL system('mkdir '//trim(filen))
      endif

      if (prlgr) then
        filen=trim(outfol)//'/IB_out'
        inquire(file=filen, exist=dir_e)
        if (.not.dir_e) then
          CALL system('mkdir '//trim(filen))
        endif
      endif

      open(141,file=trim(outfol)//'/simulation_report.out', &
              status='unknown') 

      open(115,file=trim(outfol)//'/fluid_rep.out', &
                status='unknown')


      if (ibact) then
        
        filen=trim(outfol)//'/forcecomp.out'
        open(126,file=trim(filen),status='unknown')

        filen=trim(outfol)//'/rbkin.out'
        open(39, file=filen, action="write")
        
      endif

     !write on simulation report 
      CALL date_and_time(b(1),b(2),b(3),det)

 101  format (' Start date: ',i2.2,'/',i2.2,'/',i4.4)
      write(141,101) det(3),det(2),det(1)

 102  format (' Start time: ',i2.2,':',i2.2,':',i2.2)
      write(141,102) det(5),det(6),det(7)


      ENDSUBROUTINE   
      
     !===================================================================


