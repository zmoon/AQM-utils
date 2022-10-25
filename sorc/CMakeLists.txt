# This is the cmake build file for the src directory in AQM-utils
#
add_subdirectory(gefs2clbcs_para)
add_subdirectory(decomp-pt-mpi)
if (BUILD_POST_STAT)
  add_subdirectory(O3-PM25-stat)
endif()
