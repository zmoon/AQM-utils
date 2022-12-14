help([[
This loads libraries for building NEXUS on
the NOAA operational machine WCOSS2 using Intel-19.1.3.304
]])

whatis([===[Loads libraries needed for building the UFS SRW App on WCOSS2 ]===])

load(pathJoin("envvar", os.getenv("envvar_ver")))

load(pathJoin("PrgEnv-intel", os.getenv("PrgEnv_intel_ver")))
load(pathJoin("intel", os.getenv("intel_ver")))
load(pathJoin("craype", os.getenv("craype_ver")))
load(pathJoin("cray-mpich", os.getenv("cray_mpich_ver")))

load(pathJoin("cmake", os.getenv("cmake_ver")))

setenv("HPC_OPT","/apps/ops/para/libs")
prepend_path("MODULEPATH", pathJoin("/apps/ops/para/libs/modulefiles/compiler/intel", os.getenv("intel_ver")))
prepend_path("MODULEPATH", pathJoin("/apps/ops/para/libs/modulefiles/mpi/intel", os.getenv("intel_ver"), "cray-mpich", os.getenv("cray_mpich_ver")))

load(pathJoin("cray-pals", os.getenv("cray_pals_ver")))

load(pathJoin("hdf5", os.getenv("hdf5_ver")))
load(pathJoin("netcdf", os.getenv("netcdf_ver")))
load(pathJoin("esmf", os.getenv("esmf_ver")))

setenv("CMAKE_C_COMPILER","cc")
setenv("CMAKE_CXX_COMPILER","CC")
setenv("CMAKE_Fortran_COMPILER","ftn")
setenv("CMAKE_Platform","wcoss2")
