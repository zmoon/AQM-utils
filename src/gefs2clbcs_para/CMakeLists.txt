set(exe_name gefs2lbc_para)
set(lib_name gefs2lbc_para_lib)
set(exe_src gefs2lbc_para.f90)

add_library(${lib_name} STATIC ${exe_src})
add_executable(${exe_name} ${exe_src})

set(mod_dir "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(${lib_name} PROPERTIES Fortran_MODULE_DIRECTORY ${mod_dir})
target_include_directories(${lib_name} INTERFACE ${mod_dir})

target_link_libraries(
  ${lib_name}
  PUBLIC 
  nemsio::nemsio
  bacio::bacio_4
  w3nco::w3nco_d
  MPI::MPI_Fortran
  esmf
  NetCDF::NetCDF_Fortran
)

target_link_libraries(${exe_name} PRIVATE ${lib_name})

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})
