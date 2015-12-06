#!/bin/bash

cmake \
-D CMAKE_INSTALL_PREFIX=/home/yf30/alps-2.2.b3-r7462-src-with-boost/hdf5_installed/ \
-D BUILD_SHARED_LIBS=OFF \
-D CMAKE_C_COMPILER=/usr/bin/gcc \
-D CMAKE_CXX_COMPILER=/usr/bin/g++ \
/home/yf30/HDF5-1.8.12/


#-D CMAKE_EXE_LINKER_FLAGS=-lgfortran \
