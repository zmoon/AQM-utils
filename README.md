# AQM-utils

## Build
```
cd AQM-utils
module use modulefiles
source env/build_[machine]_[compiler].env
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=..
make -j2
make install
```
