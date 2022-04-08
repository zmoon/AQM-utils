# check environment
cmake_minimum_required(VERSION 3.15)

# set compiler
if (NOT DEFINED ENV{CMAKE_Fortran_COMPILER})
  message(FATAL_ERROR "CMAKE_Fortran_COMPILER is not defined")
endif()

set(CMAKE_Fortran_COMPILER $ENV{CMAKE_Fortran_COMPILER})

# set the project name and version
project(
  AQM-utils
  VERSION 1.0.0
  LANGUAGES Fortran)

# include local macros
list (APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# find external libraries
find_package(ESMF 8.0.0 REQUIRED)
find_package(NetCDF 4.3.3 REQUIRED COMPONENTS Fortran)
find_package(MPI REQUIRED Fortran)
find_package(w3nco 2.4.0 REQUIRED)
find_package(bacio 2.4.0 REQUIRED)
find_package(nemsio 2.5.0 REQUIRED)
find_package(sigio 2.3.0 REQUIRED)
find_package(sp 2.3.3 REQUIRED)
find_package(ip 3.3.3 REQUIRED)
find_package(g2 3.4.0 REQUIRED)
find_package(sigio 2.3.0 REQUIRED)

# EMC requires executables in ./exec
set(exec_dir bin)
if(EMC_EXEC_DIR)
  set(exec_dir exec)
endif()

# add project's subdirectories
add_subdirectory(src)
