!
! === Process PM25 statistics with all the fcst hours in one file
!

      use wgrib2api

      integer :: nfcst, iret, ihrs, nhrs, ihrs1
      real :: focushrs
      real, allocatable :: grid(:,:), PDMAX24(:,:), PDMAX1(:,:)
      real, allocatable :: OZMAX1(:,:), OZMAX8(:,:), OZAVE8(:,:)
      character (len=200) :: file, inv, var, var2, var3, output
      character (len=99) :: invline, stri, svar2(72), ihrx, nhrx
      character (len=500) :: gridline

!----------------------

      file = 'PMTF_OZCON-01_2_48.grb2'   ! 'PMTF-01_2_48.grb2'
      inv = '@mem:0'                     ! stored in memory and no file

      output = 'PM_O3_stat.grb2'

! make inv file, save in memory file #0

      iret = grb2_mk_inv(file, inv)
      if (iret.ne.0) stop 1

! forcast hours for 12Z cycle to be used (04Z-04Z)

      ihrs=16
      nhrs=40

      ihrs1=ihrs+1
      focushrs=nhrs-ihrs+1

! initialize the variables

      PDMAX1 = 0.
      PDMAX24= 0
      OZMAX1 = 0.
      OZMAX8 = 0.

! Strings for searching in the "desc"

      var  = ':PMTF:'
      var2 = ' hour fcst:'

! get number of forecast hours in the file

      nfcst = grb2_inq(file,inv,var,desc=invline)

! process fcst data from the beginning

      do i = ihrs, nhrs  ! nfcst
         write(stri,*) i      ! convert number "i" to string "stri"
         svar2(i)=':'//adjustl(trim(stri)//trim(var2))
      enddo


! --- PMTF

      var  = ':PMTF:'

      do i = ihrs, nhrs
         iret = grb2_inq(file,inv,var,trim(svar2(i)),data2=grid)

         !=================================================
         !  for wgrib2, the 1st data need to be assigned 
         !  for dimension match
         !=================================================

         if (i==ihrs) then
            PDMAX24 = grid
         else
            PDMAX24 = PDMAX24 + grid
         endif

         PDMAX1 = max (grid,PDMAX1)
      enddo

      PDMAX24=PDMAX24/focushrs

! --- OZCON

      var  = ':OZCON:'

      do i = ihrs, nhrs
         iret = grb2_inq(file,inv,var,trim(svar2(i)),data2=grid)
         OZMAX1 = max (grid,OZMAX1)
      enddo

      !--- O3 8HR MAX

      do kk = 1, focushrs
         i1 = ihrs+kk-1
         i8 = i1-7               ! if forward "i1+7"

         do i= i8, i1
            iret = grb2_inq(file,inv,var,trim(svar2(i)),data2=grid)
 
            if (i==i8) then      ! if forward "i1"
               OZAVE8 = grid
            else
               OZAVE8 = OZAVE8 + grid
            endif
         enddo

         OZAVE8 = OZAVE8/8.
         OZMAX8 = max (OZAVE8,OZMAX8)
      enddo

! new method: override metadata

      write(ihrx,*) ihrs1     ! should be "ihrs", temporary used to match the MetPlus criteria
      write(nhrx,*) nhrs
      var3=adjustl(trim(ihrx)//'-'//trim(nhrx)//' hour ave fcst')

      iret = grb2_wrt(output,file,1,data2=PDMAX24,meta=invline,var='PDMAX24',timing=var3)
      if (iret.ne.0) stop 3

      iret = grb2_wrt(output,file,1,data2=PDMAX1,meta=invline,var='PDMAX1',timing=var3,append=1)
      if (iret.ne.0) stop 4

      iret = grb2_wrt(output,file,1,data2=OZMAX1,meta=invline,var='OZMAX1',timing=var3,append=1)
      if (iret.ne.0) stop 5

      iret = grb2_wrt(output,file,1,data2=OZMAX8,meta=invline,var='OZMAX8',timing=var3,append=1)
      if (iret.ne.0) stop 6


      stop
      end
