# AQM-utils

## Build
```
cd AQM-utils
module use modulefiles
module load [machine].[compiler]
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=.. -DCMAKE_INSTALL_BINDIR=exec
make -j2
make install
```
