# AQM-utils

## Build
```
cd AQM-utils
module use modulefiles
source modulefiles/build_[machine]_[compiler].env
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=.. -DCMAKE_INSTALL_BINDIR=exec
make -j2
make install
```
