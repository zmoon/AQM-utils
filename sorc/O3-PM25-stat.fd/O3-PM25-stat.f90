!
! === Process O3 and PM2.5 statistics
!
! prototyped by Hsin-mu Lin, Aug, 2022
! updated by Kai Wang and Jianping Huang, Oct, 2022

      use wgrib2api

      integer :: nfcst, iret 
      integer :: ihrs1, nhrs1, ihrs2, nhrs2, ihrs3, nhrs3
      integer :: ihrs1x, nhrs1x, ihrs2x, nhrs2x, ihrs3x, nhrs3x
      real, allocatable :: grid(:,:)
      real, allocatable :: OZMAX1(:,:), OZMAX8(:,:), OZAVE8(:,:)
      real, allocatable :: PDMAX24(:,:), PDMAX1(:,:)
      character (len=200) :: inv, var, var3
      character (len=200) :: file1_1h, file2_1h, file3_1h
      character (len=200) :: file1_8h, file2_8h, file3_8h
      character (len=200) :: file1_pm, file2_pm, file3_pm
      character (len=200) :: output1_o3, output2_o3
      character (len=200) :: output1_pm, output2_pm
      character (len=99)  :: invline, invlinex
      character (len=99)  :: ihrx1, nhrx1, ihrx2, nhrx2, ihrx3, nhrx3
      character (len=99)  :: ihrx1x, nhrx1x, ihrx2x, nhrx2x, ihrx3x, nhrx3x

!===================================================================
! ** open time step files
!===================================================================

      open(9,file='meta_time-step.txt', form='formatted')
      open(10,file='meta_time-step.max8ho3.txt', form='formatted')

!===================================================================
! ** note: input files (file1, file2, file3) are pre-processed to
!    include only data needed for daily 04z-04z statistics 
!===================================================================

      file1_1h = 'day1-POST-UPP-INPUT-OZCON.grib2'
      file2_1h = 'day2-POST-UPP-INPUT-OZCON.grib2'
      file3_1h = 'day3-POST-UPP-INPUT-OZCON.grib2'

      file1_8h = 'day1-POST-UPP-INPUT-OZCON_8hrmax.grib2'
      file2_8h = 'day2-POST-UPP-INPUT-OZCON_8hrmax.grib2'
      file3_8h = 'day3-POST-UPP-INPUT-OZCON_8hrmax.grib2'

      file1_pm = 'day1-POST-UPP-INPUT-PMTF.grib2'
      file2_pm = 'day2-POST-UPP-INPUT-PMTF.grib2'
      file3_pm = 'day3-POST-UPP-INPUT-PMTF.grib2'

      output1_o3 = 'aqm.max_1hr_o3.grib2'
      output2_o3 = 'aqm.max_8hr_o3.grib2'

      output1_pm = 'aqm.ave_24hr_pm25.grib2'
      output2_pm = 'aqm.max_1hr_pm25.grib2'

!for O3
!-----------------------

      var  = ':OZCON:'

!========================
! For day1
!========================

!=================
! === O3 daily max

! make inv file, save in memory file #0

      inv = '@mem:0'                     ! stored in memory and no file

      iret = grb2_mk_inv(file1_1h, inv)

      if (iret.ne.0) then
         ! stop 1
         go to 66  ! skip to day2 if no day1 data
      endif

! get number of forecast hours in the file, also "invline" info

      !---NOTE: "invline" is needed for the "grb2_wrt" later

      nfcst = grb2_inq(file1_1h,inv,var,desc=invline)

!=== process post data from 04Z-04Z

      OZMAX1 = 0.

      do i = 1, nfcst
         iret = grb2_inq(file1_1h,inv,var,sequential=i-1,data2=grid)
         OZMAX1 = max (grid,OZMAX1)
      enddo

!=================
!=== O3 8HR MAX

! make inv file, save in memory file #0

      inv = '@mem:1'                     ! stored in memory and no file

      iret = grb2_mk_inv(file1_8h, inv)
      if (iret.ne.0) stop 10

      !---NOTE: "invline" is needed for the "grb2_wrt" later

      nfcst = grb2_inq(file1_8h,inv,var,desc=invlinex)

      OZMAX8 = 0.

      do ii = 1, nfcst-7

         do i= 1, ii+7     ! for 8 hours
            iret = grb2_inq(file1_8h,inv,var,sequential=i-1,data2=grid)

            if (i >= ii) then
               if (i==ii) then
                  OZAVE8 = grid
               else
                  OZAVE8 = OZAVE8 + grid
               endif
            endif
         enddo

         OZMAX8 = max (OZAVE8,OZMAX8)

      enddo

      OZMAX8 = OZMAX8/8.

! write day1 metadata

      read (9,*) ihrs1
      read (9,*) nhrs1

      ! ihrs=17
      ! nhrs=16

      write(ihrx1,*) ihrs1
      write(nhrx1,*) nhrs1
      var3=adjustl(trim(ihrx1)//'-'//trim(nhrx1)//' hour ave fcst')

      iret = grb2_wrt(output1_o3,file1_1h,1,data2=OZMAX1,meta=invline,var='OZMAX1',timing=var3)
      if (iret.ne.0) stop 2

      read (10,*) ihrs1x
      read (10,*) nhrs1x

      write(ihrx1x,*) ihrs1x
      write(nhrx1x,*) nhrs1x
      var3=adjustl(trim(ihrx1x)//'-'//trim(nhrx1x)//' hour ave fcst')

      iret = grb2_wrt(output2_o3,file1_8h,1,data2=OZMAX8,meta=invlinex,var='OZMAX8',timing=var3)
      if (iret.ne.0) stop 3

66    continue

!========================
! For day2
!========================

!=================
! === O3 daily max

      inv = '@mem:2'

! make inv file, save in memory file #0

      iret = grb2_mk_inv(file2_1h, inv)
      if (iret.ne.0) stop 4

! get number of forecast hours in the file, also "invline" info

      !---NOTE: "invline" is needed for the "grb2_wrt" later

      nfcst = grb2_inq(file2_1h,inv,var,desc=invline)

!=== process post data from 04Z-04Z

      OZMAX1 = 0.

      do i = 1, nfcst
         iret = grb2_inq(file2_1h,inv,var,sequential=i-1,data2=grid)
         OZMAX1 = max (grid,OZMAX1)
      enddo

!=================
!=== O3 8HR MAX

! make inv file, save in memory file #0

      inv = '@mem:3'                     ! stored in memory and no file

      iret = grb2_mk_inv(file2_8h, inv)
      if (iret.ne.0) stop 20

      !---NOTE: "invline" is needed for the "grb2_wrt" later

      nfcst = grb2_inq(file2_8h,inv,var,desc=invlinex)

      OZMAX8 = 0.

      do ii = 1, nfcst-7

         do i= 1, ii+7     ! for 8 hours
            iret = grb2_inq(file2_8h,inv,var,sequential=i-1,data2=grid)

            if (i >= ii) then
               if (i==ii) then
                  OZAVE8 = grid
               else
                  OZAVE8 = OZAVE8 + grid
               endif
            endif
         enddo

         OZMAX8 = max (OZAVE8,OZMAX8)

      enddo

      OZMAX8 = OZMAX8/8.

! append day2 metadata

      read (9,*) ihrs2
      read (9,*) nhrs2

      write(ihrx2,*) ihrs2
      write(nhrx2,*) nhrs2
      var3=adjustl(trim(ihrx2)//'-'//trim(nhrx2)//' hour ave fcst')

      iret = grb2_wrt(output1_o3,file2_1h,1,data2=OZMAX1,meta=invline,var='OZMAX1',timing=var3,append=1)
      if (iret.ne.0) stop 5

      read (10,*) ihrs2x
      read (10,*) nhrs2x

      write(ihrx2x,*) ihrs2x
      write(nhrx2x,*) nhrs2x
      var3=adjustl(trim(ihrx2x)//'-'//trim(nhrx2x)//'hour ave fcst')

      iret = grb2_wrt(output2_o3,file2_8h,1,data2=OZMAX8,meta=invlinex,var='OZMAX8',timing=var3,append=1)
      if (iret.ne.0) stop 6

!========================
! For day3
!========================

!=================
! === O3 daily max

      inv = '@mem:4'

! make inv file, save in memory file #0

      iret = grb2_mk_inv(file3_1h, inv)

      if (iret.ne.0) then
         ! stop 7
         go to 77 ! skip to end if no day3 data
      endif

! get number of forecast hours in the file, also "invline" info

      !---NOTE: "invline" is needed for the "grb2_wrt" later

      nfcst = grb2_inq(file3_1h,inv,var,desc=invline)

!=== process post data from 04Z-04Z

      OZMAX1 = 0.

      do i = 1, nfcst
         iret = grb2_inq(file3_1h,inv,var,sequential=i-1,data2=grid)
         OZMAX1 = max (grid,OZMAX1)
      enddo

!=================
!=== O3 8HR MAX

! make inv file, save in memory file #0

      inv = '@mem:5'                     ! stored in memory and no file

      iret = grb2_mk_inv(file3_8h, inv)
      if (iret.ne.0) stop 30

      !---NOTE: "invline" is needed for the "grb2_wrt" later

      nfcst = grb2_inq(file3_8h,inv,var,desc=invlinex)

      OZMAX8 = 0.

      do ii = 1, nfcst-7

         do i= 1, ii+7     ! for 8 hours
            iret = grb2_inq(file3_8h,inv,var,sequential=i-1,data2=grid)

            if (i >= ii) then
               if (i==ii) then
                  OZAVE8 = grid
               else
                  OZAVE8 = OZAVE8 + grid
               endif
            endif
         enddo

         OZMAX8 = max (OZAVE8,OZMAX8)

      enddo

      OZMAX8 = OZMAX8/8.

! append day3 metadata

      read (9,*) ihrs3
      read (9,*) nhrs3

      ! ihrs=41
      ! nhrs=64

      write(ihrx3,*) ihrs3     ! should be "ihrs", temporary used toatch the MetPlus criteria
      write(nhrx3,*) nhrs3
      var3=adjustl(trim(ihrx3)//'-'//trim(nhrx3)//' hour ave fcst')

      iret = grb2_wrt(output1_o3,file3_1h,1,data2=OZMAX1,meta=invline,var='OZMAX1',timing=var3,append=1)
      if (iret.ne.0) stop 8

      read (10,*) ihrs3x
      read (10,*) nhrs3x

      write(ihrx3x,*) ihrs3x     ! should be "ihrs", temporary used toatch the MetPlus criteria
      write(nhrx3x,*) nhrs3x
      var3=adjustl(trim(ihrx3x)//'-'//trim(nhrx3x)//' hour ave fcst')

      iret = grb2_wrt(output2_o3,file3_8h,1,data2=OZMAX8,meta=invlinex,var='OZMAX8',timing=var3,append=1)
      if (iret.ne.0) stop 9

77    continue


! for PM2.5
!-----------------------

      var  = ':PMTF:'

!========================
! For day1
!========================

! make inv file, save in memory file #0

      inv = '@mem:6'                     ! stored in memory and no file

      iret = grb2_mk_inv(file1_pm, inv)

      if (iret.ne.0) then
         ! stop 1
         go to 88  ! skip to day2 if no day1 data
      endif

! get number of forecast hours in the file, also "invline" info

      !---NOTE: "invline" is needed for the "grb2_wrt" later

      nfcst = grb2_inq(file1_pm,inv,var,desc=invline)

!=== process post data from 04Z-04Z

! --- PMTF

      PDMAX1 = 0.
      PDMAX24= 0.

      do i = 1, nfcst
         iret = grb2_inq(file1_pm,inv,var,sequential=i-1,data2=grid)

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

      !read (9,*) ihrs
      !read (9,*) nhrs

      ! ihrs=17
      ! nhrs=16

      !write(ihrx,*) ihrs
      !write(nhrx,*) nhrs
      var3=adjustl(trim(ihrx1)//'-'//trim(nhrx1)//' hour ave fcst')

      iret = grb2_wrt(output1_pm,file1_pm,1,data2=PDMAX24,meta=invline,var='PMTF',timing=var3)
      if (iret.ne.0) stop 12

      iret = grb2_wrt(output2_pm,file1_pm,1,data2=PDMAX1,meta=invline,var='PDMAX1',timing=var3)
      if (iret.ne.0) stop 13

88    continue

!========================
! For day2
!========================

      inv = '@mem:7'

! make inv file, save in memory file #0

      iret = grb2_mk_inv(file2_pm, inv)
      if (iret.ne.0) stop 14

! get number of forecast hours in the file, also "invline" info

      !---NOTE: "invline" is needed for the "grb2_wrt" later

      nfcst = grb2_inq(file2_pm,inv,var,desc=invline)

!=== process post data from 04Z-04Z

      PDMAX1 = 0.
      PDMAX24= 0

      do i = 1, nfcst
         iret = grb2_inq(file2_pm,inv,var,sequential=i-1,data2=grid)

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

      !read (9,*) ihrs
      !read (9,*) nhrs

      ! ihrs=17
      ! nhrs=40

      !write(ihrx,*) ihrs
      !write(nhrx,*) nhrs
      var3=adjustl(trim(ihrx2)//'-'//trim(nhrx2)//' hour ave fcst')

      iret = grb2_wrt(output1_pm,file2_pm,1,data2=PDMAX24,meta=invline,var='PMTF',timing=var3,append=1)
      if (iret.ne.0) stop 15

      iret = grb2_wrt(output2_pm,file2_pm,1,data2=PDMAX1,meta=invline,var='PDMAX1',timing=var3,append=1)
      if (iret.ne.0) stop 16

!========================
! For day3
!========================

      inv = '@mem:8'

! make inv file, save in memory file #0

      iret = grb2_mk_inv(file3_pm, inv)

      if (iret.ne.0) then
         ! stop 7
         go to 99 ! skip to end if no day3 data
      endif

! get number of forecast hours in the file, also "invline" info

      !---NOTE: "invline" is needed for the "grb2_wrt" later

      nfcst = grb2_inq(file3_pm,inv,var,desc=invline)

!=== process post data from 04Z-04Z

      PDMAX1 = 0.
      PDMAX24= 0

      do i = 1, nfcst
         iret = grb2_inq(file3_pm,inv,var,sequential=i-1,data2=grid)

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

      !read (9,*) ihrs
      !read (9,*) nhrs

      ! ihrs=41
      ! nhrs=64

      !write(ihrx,*) ihrs
      !write(nhrx,*) nhrs
      var3=adjustl(trim(ihrx3)//'-'//trim(nhrx3)//' hour ave fcst')

      iret = grb2_wrt(output1_pm,file3_pm,1,data2=PDMAX24,meta=invline,var='PMTF',timing=var3,append=1)
      if (iret.ne.0) stop 18

      iret = grb2_wrt(output2_pm,file3_pm,1,data2=PDMAX1,meta=invline,var='PDMAX1',timing=var3,append=1)
      if (iret.ne.0) stop 19

99    continue

      close(9)
      close(10) 

      stop
      end
