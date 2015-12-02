#!/bin/bash

if [ "$1" == "static" ]; then

  echo -e "cmake: build a static ALPS library..."
  cmake \
  -D CMAKE_EXE_LINKER_FLAGS="${CMAKE_EXE_LINKER_FLAGS} -lgfortran" \
  -D CMAKE_FIND_LIBRARY_SUFFIXES=".a" \
  -D ALPS_BUILD_TESTS=OFF \
  -D ALPS_BUILD_EXAMPLES=OFF \
  -D ALPS_BUILD_APPLICATIONS=OFF \
  -D ALPS_BUILD_PYTHON=OFF \
  -D BLAS_LIBRARY=/usr/lib64/libblas.a \
  -D LAPACK_LIBRARY=/usr/lib64/liblapack.a \
  -D HDF5_LIBRARIES=/home/yf30/alps-2.2.b3-r7462-src-with-boost/hdf5_installed/lib/libhdf5.a \
  -D CMAKE_INSTALL_PREFIX=/home/yf30/alps_static/ \
  -D LAPACK_64_BIT=ON \
  -D HDF5_INCLUDE_DIR=/home/yf30/alps-2.2.b3-r7462-src-with-boost/hdf5_installed/include/ \
  -D Boost_ROOT_DIR=/home/yf30/alps-2.2.b3-r7462-src-with-boost/boost/ \
  -D BUILD_SHARED_LIBS=OFF \
  -D PYTHON_LIBRARY=/usr/lib64/libpython2.6.so \
  -D CMAKE_C_COMPILER=/usr/bin/gcc \
  -D CMAKE_CXX_COMPILER=/usr/bin/g++ \
  -D ALPS_ENABLE_MPI=ON \
  -D MPI_CXX_COMPILER=/usr/lib64/openmpi/bin/mpicxx \
  -D MPI_C_INCLUDE_PATH=/usr/include/openmpi-x86_64/ \
  -D MPI_C_LIBRARIES="/usr/lib64/openmpi/lib/libmpi.so;/usr/lib64/openmpi/lib/libopen-pal.so;/usr/lib64/openmpi/lib/libopen-rte.so" \
  -D MPI_CXX_INCLUDE_PATH=/usr/include/openmpi-x86_64/ \
  -D MPI_CXX_LIBRARIES="/usr/lib64/openmpi/lib/libmpi_cxx.so;/usr/lib64/openmpi/lib/libopen-pal.so;/usr/lib64/openmpi/lib/libopen-rte.so" \
  /home/yf30/alps-2.2.b3-r7462-src-with-boost/alps/

elif [ "$1" == "shared" ]; then

  echo -e "cmake: build a shared (dynamical) ALPS library..."
  cmake \
  -D ALPS_BUILD_TESTS=ON \
  -D ALPS_BUILD_EXAMPLES=ON \
  -D ALPS_BUILD_APPLICATIONS=ON \
  -D ALPS_BUILD_PYTHON=ON \
  -D BLAS_LIBRARY=/usr/lib64/libblas.so \
  -D LAPACK_LIBRARY=/usr/lib64/liblapack.so \
  -D CMAKE_INSTALL_PREFIX=/home/yf30/alps/ \
  -D LAPACK_64_BIT=ON \
  -D Boost_ROOT_DIR=/home/yf30/alps-2.2.b3-r7462-src-with-boost/boost/ \
  -D HDF5_LIBRARIES=/usr/lib64/libhdf5.so.6 \
  -D HDF5_INCLUDE_DIR=/home/yf30/alps-2.2.b3-r7462-src-with-boost/hdf5_installed/include/ \
  -D PYTHON_LIBRARY=/usr/lib64/libpython2.6.so \
  -D CMAKE_C_COMPILER=/usr/bin/gcc \
  -D CMAKE_CXX_COMPILER=/usr/bin/g++ \
  -D ALPS_ENABLE_MPI=ON \
  -D MPI_CXX_COMPILER=/usr/lib64/openmpi/bin/mpicxx \
  -D MPI_C_INCLUDE_PATH=/usr/include/openmpi-x86_64/ \
  -D MPI_C_LIBRARIES="/usr/lib64/openmpi/lib/libmpi.so;/usr/lib64/openmpi/lib/libopen-pal.so;/usr/lib64/openmpi/lib/libopen-rte.so" \
  -D MPI_CXX_INCLUDE_PATH=/usr/include/openmpi-x86_64/ \
  -D MPI_CXX_LIBRARIES="/usr/lib64/openmpi/lib/libmpi_cxx.so;/usr/lib64/openmpi/lib/libopen-pal.so;/usr/lib64/openmpi/lib/libopen-rte.so" \
  /home/yf30/alps-2.2.b3-r7462-src-with-boost/alps/

else

  echo -e "To build the ALPS library, this script takes one argument (\"static\" or \"shared\")."

fi


#-D CMAKE_EXE_LINKER_FLAGS="${CMAKE_EXE_LINKER_FLAGS} -lgfortran" \
#-D CMAKE_EXE_LINKER_FLAGS=-lgfortran \
