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

# add project's subdirectories
add_subdirectory(sorc)
