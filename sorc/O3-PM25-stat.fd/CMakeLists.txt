set(exe_name O3-PM25-stat)
set(lib_name O3-PM25-stat_lib)
set(exe_src O3-PM25-stat.f90)
set(CMAKE_Fortran_FLAGS "-no-wrap-margin -qopenmp")

add_library(${lib_name} STATIC ${exe_src})
add_executable(${exe_name} ${exe_src})

set(mod_dir "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(${lib_name} PROPERTIES Fortran_MODULE_DIRECTORY ${mod_dir})
target_include_directories(${lib_name} INTERFACE ${mod_dir})

# find external libraries
find_package(wgrib2 REQUIRED)

target_link_libraries(
  ${lib_name}
  PUBLIC 
  wgrib2::wgrib2
)

target_link_libraries(${exe_name} PRIVATE ${lib_name})

install(TARGETS ${exe_name})
