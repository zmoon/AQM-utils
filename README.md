# AQM-utils

## Build

- WCOSS2 (Cactus/Dogwood)
```
cd AQM-utils
source versions/build.ver.wcoss2
module use modulefiles
module load build_wcoss2.intel.lua
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=.. -DCMAKE_INSTALL_BINDIR=exec -DBUILD_POST_STAT=on
make -j2
make install
```

- Hera/Orion
```
cd AQM-utils
module use modulefiles
module load build_[machine].intel.lua
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=.. -DCMAKE_INSTALL_BINDIR=exec
make -j2
make install
```
