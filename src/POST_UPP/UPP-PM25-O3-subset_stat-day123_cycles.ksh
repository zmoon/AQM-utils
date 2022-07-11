#!/bin/ksh
#

set -xv

# ---- Get grib data at certain record # ----------------------
#
#  subset $field & ":1 hybrid level" and append all time together
#

#=======================================================================
# Directory of source script
#=======================================================================

 SDIR=$PWD

#-----------------------------------------------------------------
# get date & Cycle
#-----------------------------------------------------------------

PDY=`cut -c 1-10 cycledate`

CYCLES=`echo $PDY | cut -c 1-10`

export CYC=`echo $PDY | cut -c 9-10`

# CYCLESx=`$NDATE -24 $CYCLES`

CYCLESx24=`ndate -24 $CYCLES`
CYCLESx18=`ndate -18 $CYCLES`
CYCLESx12=`ndate -12 $CYCLES`
CYCLESx06=`ndate -06 $CYCLES`

#==================================================================

name=POST-UPP-INPUT

export field1=PMTF
export field2=OZCON
export field3='1 hybrid level'

export fname1=day1-${name}-${field1}.grib2
export fname1s=day1-${name}-${field2}.grib2
export fname1x=day1-${name}-${field2}_8hrmax.grib2

export fname2=day2-${name}-${field1}.grib2
export fname2s=day2-${name}-${field2}.grib2
export fname2x=day2-${name}-${field2}_8hrmax.grib2

export fname3=day3-${name}-${field1}.grib2
export fname3s=day3-${name}-${field2}.grib2
export fname3x=day3-${name}-${field2}_8hrmax.grib2

#----------------------------------------------------------
#---- locate the post data and get the necessary 04Z-04Z
#----------------------------------------------------------

postdata=/lfs/h2/emc/ptmp/jianping.huang/para/com/aqm/v7.0/aqm.v7.0.a1

export postdir=${postdata}/${CYCLES}/postprd

export postdirx24=${postdata}/${CYCLESx24}/postprd
export postdirx18=${postdata}/${CYCLESx18}/postprd
export postdirx12=${postdata}/${CYCLESx12}/postprd
export postdirx06=${postdata}/${CYCLESx06}/postprd

#========================================

if [ $CYC -eq 18 ]; then
  if [ -d $postdirx24 ]; then
    if [ -d $postdirx18 ]; then
      if [ -d $postdirx12 ]; then
        if [ -d $postdirx06 ]; then
          UPP-PM25-O3-subset_stat-${CYC}Z.ksh
	fi
      fi
    fi
  else
    UPP-PM25-O3-subset_stat-${CYC}Z_2_${CYC}Z.ksh
  fi
fi

if [ $CYC -eq 12 ]; then
  if [ -d $postdirx18 ]; then
    if [ -d $postdirx12 ]; then
      if [ -d $postdirx06 ]; then
        UPP-PM25-O3-subset_stat-${CYC}Z.ksh
      fi
    fi
  else  
    if [ -d $postdirx24 ]; then
      UPP-PM25-O3-subset_stat-${CYC}Z_2_${CYC}Z.ksh
    fi
  fi 
fi

if [ $CYC -eq 6 ]; then
  if [ -d $postdirx12 ]; then
    if [ -d $postdirx06 ]; then
      UPP-PM25-O3-subset_stat-${CYC}Z.ksh
    fi
  else
    if [ -d $postdirx24 ]; then
      UPP-PM25-O3-subset_stat-${CYC}Z_2_${CYC}Z.ksh
    fi
  fi
fi

if [ $CYC -eq 00 ]; then
  if [ -d $postdirx06 ]; then
    UPP-PM25-O3-subset_stat-${CYC}Z.ksh
  else
    if [ -d $postdirx24 ]; then
      UPP-PM25-O3-subset_stat-${CYC}Z_2_${CYC}Z.kshxxx
    fi
  fi
fi

#--------------------------------------

#################################################
## ---- excute the stat -----

 ./PM25-stat
 ./O3-stat

exit
