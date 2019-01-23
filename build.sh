#!/bin/bash
cd glpk/glpk-4.47
./configure
make
cd ../../ctbn_build_dir
cmake ..
# use -jx to run x compiles at the same time (for multiple processors)
# eg: replace below with "make -j4"
make
make test
cd ..
