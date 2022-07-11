!
! === Process O3 statistics
!

      use wgrib2api

      integer :: nfcst, iret, ihrs, nhrs
      real, allocatable :: grid(:,:)
      real, allocatable :: OZMAX1(:,:), OZMAX8(:,:), OZAVE8(:,:)
      character (len=200) :: inv, var, var3
      character (len=200) :: file1, file2, file3
      character (len=200) :: file1x, file2x, file3x
      character (len=200) :: output1, output2
      character (len=99)  :: invline, ihrx, nhrx, invlinex

!===================================================================
! ** open time step file
!===================================================================

      open(9,file='meta_time-step.txt', form='formatted')

!===================================================================
! ** note: input files (file1, file2, file3) are pre-processed to
!    include only data needed for daily 04z-04z statistics 
!===================================================================

      file1 = 'day1-POST-UPP-INPUT-OZCON.grib2'
      file2 = 'day2-POST-UPP-INPUT-OZCON.grib2'
      file3 = 'day3-POST-UPP-INPUT-OZCON.grib2'

      file1x= 'day1-POST-UPP-INPUT-OZCON_8hrmax.grib2'
      file2x= 'day2-POST-UPP-INPUT-OZCON_8hrmax.grib2'
      file3x= 'day3-POST-UPP-INPUT-OZCON_8hrmax.grib2'

      output1 = 'aqm.t12z.max_1hr_o3.grib2'
      output2 = 'aqm.t12z.max_8hr_o3.grib2'

!-----------------------

      var  = ':OZCON:'

!========================
! For day1
!========================

!=================
! === O3 daily max

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

      OZMAX1 = 0.

      do i = 1, nfcst
         iret = grb2_inq(file1,inv,var,sequential=i-1,data2=grid)
         OZMAX1 = max (grid,OZMAX1)
      enddo

!=================
!=== O3 8HR MAX

! make inv file, save in memory file #0

      inv = '@mem:1'                     ! stored in memory and no file

      iret = grb2_mk_inv(file1x, inv)
      if (iret.ne.0) stop 10

      !---NOTE: "invline" is needed for the "grb2_wrt" later

      nfcst = grb2_inq(file1x,inv,var,desc=invlinex)

      OZMAX8 = 0.

      do ii = 1, nfcst-7

         do i= 1, ii+7     ! for 8 hours
            iret = grb2_inq(file1x,inv,var,sequential=i-1,data2=grid)
            !if (iret==0) go to 66

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

66    OZMAX8 = OZMAX8/8.

! write day1 metadata

      read (9,*) ihrs
      read (9,*) nhrs

      ! ihrs=17
      ! nhrs=16

      write(ihrx,*) ihrs
      write(nhrx,*) nhrs
      var3=adjustl(trim(ihrx)//'-'//trim(nhrx)//' hour ave fcst')

      iret = grb2_wrt(output1,file1,1,data2=OZMAX1,meta=invline,var='OZMAX1',timing=var3)
      if (iret.ne.0) stop 2

      iret = grb2_wrt(output2,file1x,1,data2=OZMAX8,meta=invlinex,var='OZMAX8',timing=var3)
      if (iret.ne.0) stop 3

77    continue

!========================
! For day2
!========================

!=================
! === O3 daily max

      inv = '@mem:2'

! make inv file, save in memory file #0

      iret = grb2_mk_inv(file2, inv)
      if (iret.ne.0) stop 4

! get number of forecast hours in the file, also "invline" info

      !---NOTE: "invline" is needed for the "grb2_wrt" later

      nfcst = grb2_inq(file2,inv,var,desc=invline)

!=== process post data from 04Z-04Z

      OZMAX1 = 0.

      do i = 1, nfcst
         iret = grb2_inq(file2,inv,var,sequential=i-1,data2=grid)
         OZMAX1 = max (grid,OZMAX1)
      enddo

!=================
!=== O3 8HR MAX

! make inv file, save in memory file #0

      inv = '@mem:3'                     ! stored in memory and no file

      iret = grb2_mk_inv(file2x, inv)
      if (iret.ne.0) stop 20

      !---NOTE: "invline" is needed for the "grb2_wrt" later

      nfcst = grb2_inq(file2x,inv,var,desc=invlinex)

      OZMAX8 = 0.

      do ii = 1, nfcst-7

         do i= 1, ii+7     ! for 8 hours
            iret = grb2_inq(file2x,inv,var,sequential=i-1,data2=grid)
            ! if (iret==0) go to 33

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

33    OZMAX8 = OZMAX8/8.

! append day2 metadata

      read (9,*) ihrs
      read (9,*) nhrs

      ! ihrs=17
      ! nhrs=40

      write(ihrx,*) ihrs
      write(nhrx,*) nhrs
      var3=adjustl(trim(ihrx)//'-'//trim(nhrx)//' hour ave fcst')

      iret = grb2_wrt(output1,file2,1,data2=OZMAX1,meta=invline,var='OZMAX1',timing=var3,append=1)
      if (iret.ne.0) stop 5

      iret = grb2_wrt(output2,file2x,1,data2=OZMAX8,meta=invlinex,var='OZMAX8',timing=var3,append=1)
      if (iret.ne.0) stop 6

!========================
! For day3
!========================

!=================
! === O3 daily max

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

      OZMAX1 = 0.

      do i = 1, nfcst
         iret = grb2_inq(file3,inv,var,sequential=i-1,data2=grid)
         OZMAX1 = max (grid,OZMAX1)
      enddo

!=================
!=== O3 8HR MAX

! make inv file, save in memory file #0

      inv = '@mem:5'                     ! stored in memory and no file

      iret = grb2_mk_inv(file3x, inv)
      if (iret.ne.0) stop 30

      !---NOTE: "invline" is needed for the "grb2_wrt" later

      nfcst = grb2_inq(file3x,inv,var,desc=invlinex)

      OZMAX8 = 0.

      do ii = 1, nfcst-7

         do i= 1, ii+7     ! for 8 hours
            iret = grb2_inq(file3x,inv,var,sequential=i-1,data2=grid)
            ! if (iret==0) go to 22

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

22    OZMAX8 = OZMAX8/8.

! append day3 metadata

      read (9,*) ihrs
      read (9,*) nhrs

      ! ihrs=41
      ! nhrs=64

      write(ihrx,*) ihrs     ! should be "ihrs", temporary used toatch the MetPlus criteria
      write(nhrx,*) nhrs
      var3=adjustl(trim(ihrx)//'-'//trim(nhrx)//' hour ave fcst')

      iret = grb2_wrt(output1,file3,1,data2=OZMAX1,meta=invline,var='OZMAX1',timing=var3,append=1)
      if (iret.ne.0) stop 8

      iret = grb2_wrt(output2,file3x,1,data2=OZMAX8,meta=invlinex,var='OZMAX8',timing=var3,append=1)
      if (iret.ne.0) stop 9

88    continue

      close(9) 

      stop
      end
