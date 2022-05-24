# General Bspline Representation Library (libbspline)

## Installation
libbspline is a header only library, several test nodes are compiled using CMake.

### Setup
```bash
cd <libbspline directory>/cpp
mkdir build && cd build
cmake .. 
make
```

#### CTest
To run sample scripts to synchronize MATLAB and C++ output
```bash
make test
```

Where the CTest output in the console should be as follows
```bash
Running tests...
Test project /home/nvidia/libbspline/cpp/build
    Start 1: m_matrix
1/4 Test #1: m_matrix .........................   Passed    0.00 sec
    Start 2: txt_reader
2/4 Test #2: txt_reader .......................   Passed    0.00 sec
    Start 3: 1d_bspline
3/4 Test #3: 1d_bspline .......................   Passed    0.00 sec
    Start 4: 3d_bspline
4/4 Test #4: 3d_bspline .......................   Passed    0.00 sec

100% tests passed, 0 tests failed out of 4

Total Test time (real) =   0.01 sec
```

#### Include in other projects:
To link this lib properly, add following in the `CMakeLists.txt`
```
find_package(libbspline REQUIRED)
include_directories(${LIBBSPLINE_INCLUDE_DIRS})
```

### References
1. Liu Sikang's Decomp Util https://github.com/sikang/DecompUtil