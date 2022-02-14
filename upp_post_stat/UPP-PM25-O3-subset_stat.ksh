#!/bin/ksh
#

set -xv

# ---- Get grib data at certain record # ----------------------
#
#  subset $field & ":1 hybrid level" and append all time together
#

let istart=1
let iend=48

typeset -Z2 istart
typeset -Z2 iend

name=${istart}_2_${iend}.grb2

field1=PMTF
field2=OZCON
field3='1 hybrid level'

fname=${field1}_${field2}

CYCLES=2019070912

let i=$istart

#---- location of post data

postdir=/gpfs/dell2/ptmp/Hsin-Mu.Lin/zzz-13km-new/post_fv3_${CYCLES}

while [ $i -le $iend ]; do
  let i2=i
  typeset -Z2 i2
  file=NATLEV${i2}.tm00

   if [ $i2 -eq $istart ]; then
     wgrib2 $postdir/$file -match ":$field1" -match ":$field3" -GRIB  $fname-$name
     wgrib2 $postdir/$file -match ":$field2" -match ":$field3" -append -GRIB $fname-$name
   else
     wgrib2 $postdir/$file -match ":$field1" -match ":$field3" -append -GRIB $fname-$name
     wgrib2 $postdir/$file -match ":$field2" -match ":$field3" -append -GRIB $fname-$name
   fi

  let i=i+1
done      # while

## ---- excute the stat -----

  /gpfs/dell2/emc/modeling/save/Hsin-Mu.Lin/RRFS-CMAQ-UPP_stat/PM25-O3-stat

exit
