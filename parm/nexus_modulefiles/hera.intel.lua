help([[
This module loads libraries for building NEXUS on
the NOAA RDHPC machine Hera using Intel-2022.1.2
]])

whatis([===[Loads libraries needed for building the UFS SRW App on Hera ]===])

prepend_path("MODULEPATH","/contrib/sutils/modulefiles")
load("sutils")

load(pathJoin("cmake", os.getenv("cmake_ver") or "3.20.1"))

intel_ver=os.getenv("intel_ver") or "2022.1.2"
load(pathJoin("intel", intel_ver))

impi_ver=os.getenv("impi_ver") or "2022.1.2"
load(pathJoin("impi", impi_ver))

prepend_path("MODULEPATH","/scratch1/NCEPDEV/nems/role.epic/hpc-stack/libs/intel-2022.1.2/modulefiles/stack")

load(pathJoin("hpc", os.getenv("hpc_ver") or "1.2.0"))
load(pathJoin("hpc-intel", os.getenv("hpc_intel_ver") or "2022.1.2"))
load(pathJoin("hpc-impi", os.getenv("hpc_impi_ver") or "2022.1.2"))

load(pathJoin("hdf5", os.getenv("hdf5_ver") or "1.10.6"))
load(pathJoin("netcdf", os.getenv("netcdf_ver") or "4.7.4"))
load(pathJoin("esmf", os.getenv("esmf_ver") or "8.2.1b04"))

setenv("CMAKE_C_COMPILER","mpiicc")
setenv("CMAKE_CXX_COMPILER","mpiicpc")
setenv("CMAKE_Fortran_COMPILER","mpiifort")
setenv("CMAKE_Platform","hera.intel")

