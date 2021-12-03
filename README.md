# f4d

A self-contained version of FILO for DIMACS.

# Build
```
git clone git@github.com:falcopt/f4d.git
cd f4d
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_DIMACS=1
make -j
```