!
! === Process PM25 statistics
!

      use wgrib2api

      integer :: nfcst, iret, ihrs, nhrs
      real, allocatable :: grid(:,:), PDMAX24(:,:), PDMAX1(:,:)
      character (len=200) :: inv, var, var3
      character (len=200) :: file1, file2, file3
      character (len=200) :: output1, output2
      character (len=99)  :: invline, ihrx, nhrx

!===================================================================
! ** open time step file
!===================================================================

      open(9,file='meta_time-step.txt', form='formatted')

!===========================================================================
! ** note: input files (file1, file2, file3) are pre-processed to include
!    only 04z-04z data for each day
!===========================================================================

      file1 = 'day1-POST-UPP-INPUT-PMTF.grib2'
      file2 = 'day2-POST-UPP-INPUT-PMTF.grib2'
      file3 = 'day3-POST-UPP-INPUT-PMTF.grib2'

      output1 = 'aqm.t12z.ave_24hr_pm25.grib2'
      output2 = 'aqm.t12z.max_1hr_pm25.grib2'

!-----------------------

      var  = ':PMTF:'

!========================
! For day1
!========================

! make inv file, save in memory file #0

      inv = '@mem:0'                     ! stored in memory and no file

      iret = grb2_mk_inv(file1, inv)

      if (iret.ne.0) then
         ! stop 1
         go to 77  ! skip to day2 if no day1 data
      endif

! get number of forecast hours in the file, also "invline" info

      !---NOTE: "invline" is needed for the "grb2_wrt" later

      nfcst = grb2_inq(file1,inv,var,desc=invline)

!=== process post data from 04Z-04Z

! --- PMTF

      PDMAX1 = 0.
      PDMAX24= 0.

      do i = 1, nfcst
         iret = grb2_inq(file1,inv,var,sequential=i-1,data2=grid)

         !=================================================
         !  for wgrib2, the 1st data need to be assigned
         !  for dimension match
         !=================================================

         if (i==1) then
            PDMAX24 = grid
         else
            PDMAX24 = PDMAX24 + grid
         endif

         PDMAX1 = max (grid,PDMAX1)
      enddo

      PDMAX24=PDMAX24/nfcst

! write day1 metadata

      read (9,*) ihrs
      read (9,*) nhrs

      ! ihrs=17
      ! nhrs=16

      write(ihrx,*) ihrs
      write(nhrx,*) nhrs
      var3=adjustl(trim(ihrx)//'-'//trim(nhrx)//' hour ave fcst')

      iret = grb2_wrt(output1,file1,1,data2=PDMAX24,meta=invline,var='PDMAX24',timing=var3)
      if (iret.ne.0) stop 2

      iret = grb2_wrt(output2,file1,1,data2=PDMAX1,meta=invline,var='PDMAX1',timing=var3)
      if (iret.ne.0) stop 3

77    continue

!========================
! For day2
!========================

      inv = '@mem:2'

! make inv file, save in memory file #0

      iret = grb2_mk_inv(file2, inv)
      if (iret.ne.0) stop 4

! get number of forecast hours in the file, also "invline" info

      !---NOTE: "invline" is needed for the "grb2_wrt" later

      nfcst = grb2_inq(file2,inv,var,desc=invline)

!=== process post data from 04Z-04Z

      PDMAX1 = 0.
      PDMAX24= 0

      do i = 1, nfcst
         iret = grb2_inq(file2,inv,var,sequential=i-1,data2=grid)

         !=================================================
         !  for wgrib2, the 1st data need to be assigned 
         !  for dimension match
         !=================================================

         if (i==1) then
            PDMAX24 = grid
         else
            PDMAX24 = PDMAX24 + grid
         endif

         PDMAX1 = max (grid,PDMAX1)
      enddo

      PDMAX24=PDMAX24/nfcst

! append day2 metadata

      read (9,*) ihrs
      read (9,*) nhrs

      ! ihrs=17
      ! nhrs=40

      write(ihrx,*) ihrs
      write(nhrx,*) nhrs
      var3=adjustl(trim(ihrx)//'-'//trim(nhrx)//' hour ave fcst')

      iret = grb2_wrt(output1,file2,1,data2=PDMAX24,meta=invline,var='PDMAX24',timing=var3,append=1)
      if (iret.ne.0) stop 5

      iret = grb2_wrt(output2,file2,1,data2=PDMAX1,meta=invline,var='PDMAX1',timing=var3,append=1)
      if (iret.ne.0) stop 6

!========================
! For day3
!========================

      inv = '@mem:4'

! make inv file, save in memory file #0

      iret = grb2_mk_inv(file3, inv)

      if (iret.ne.0) then
         ! stop 7
         go to 88 ! skip to end if no day3 data
      endif

! get number of forecast hours in the file, also "invline" info

      !---NOTE: "invline" is needed for the "grb2_wrt" later

      nfcst = grb2_inq(file3,inv,var,desc=invline)

!=== process post data from 04Z-04Z

      PDMAX1 = 0.
      PDMAX24= 0

      do i = 1, nfcst
         iret = grb2_inq(file3,inv,var,sequential=i-1,data2=grid)

         !=================================================
         !  for wgrib2, the 1st data need to be assigned
         !  for dimension match
         !=================================================

         if (i==1) then
            PDMAX24 = grid
         else
            PDMAX24 = PDMAX24 + grid
         endif

         PDMAX1 = max (grid,PDMAX1)
      enddo

      PDMAX24=PDMAX24/nfcst

! append day3 metadata

      read (9,*) ihrs
      read (9,*) nhrs

      ! ihrs=41
      ! nhrs=64

      write(ihrx,*) ihrs
      write(nhrx,*) nhrs
      var3=adjustl(trim(ihrx)//'-'//trim(nhrx)//' hour ave fcst')

      iret = grb2_wrt(output1,file3,1,data2=PDMAX24,meta=invline,var='PDMAX24',timing=var3,append=1)
      if (iret.ne.0) stop 8

      iret = grb2_wrt(output2,file3,1,data2=PDMAX1,meta=invline,var='PDMAX1',timing=var3,append=1)
      if (iret.ne.0) stop 9

88    continue

      close(9)
      stop
      end
