#!/bin/bash

export LIBRARY_PATH=$LIBRARY_PATH:/usr/lib64/openmpi/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64/openmpi/lib/

if [ "$1" == "static" ]; then

  # assumed static HDF5 (libhdf5.a) is installed
  # without the flag -lgfortran there will be errors (don't know why)
  echo -e "cmake: build a static ALPS library..."
  cmake \
  -D CMAKE_FIND_LIBRARY_SUFFIXES=".a" \
  -D ALPS_BUILD_TESTS=OFF \
  -D ALPS_BUILD_EXAMPLES=OFF \
  -D ALPS_BUILD_APPLICATIONS=OFF \
  -D ALPS_BUILD_PYTHON=OFF \
  -D CMAKE_INSTALL_PREFIX=/home/yf30/alps_static/ \
  -D LAPACK_64_BIT=ON \
  -D Boost_ROOT_DIR=/home/yf30/alps-2.2.b3-r7462-src-with-boost/boost/ \
  -D BUILD_SHARED_LIBS=OFF \
  -D CMAKE_C_COMPILER=/usr/bin/gcc \
  -D CMAKE_CXX_COMPILER=/usr/bin/g++ \
  -D ALPS_ENABLE_MPI=ON \
  -D BLAS_LIBRARY=/usr/lib64/libblas.a \
  -D LAPACK_LIBRARY="/usr/lib64/liblapack.a;/usr/lib/gcc/x86_64-redhat-linux/4.4.4/libgfortran.a" \
  -D HDF5_LIBRARIES="/usr/lib64/libhdf5.a;/usr/lib64/libhdf5_hl.a" \
  /home/yf30/alps-2.2.b3-r7462-src-with-boost/alps/

##  Working version (on lqcd28); don't touch this part!
##
##  assumed static HDF5 (libhdf5.a) is NOT installed
#  cmake \
#  -D CMAKE_EXE_LINKER_FLAGS="${CMAKE_EXE_LINKER_FLAGS} -lgfortran" \
#  -D CMAKE_FIND_LIBRARY_SUFFIXES=".a" \
#  -D ALPS_BUILD_TESTS=OFF \
#  -D ALPS_BUILD_EXAMPLES=OFF \
#  -D ALPS_BUILD_APPLICATIONS=OFF \
#  -D ALPS_BUILD_PYTHON=OFF \
#  -D BLAS_LIBRARY=/usr/lib64/libblas.a \
#  -D LAPACK_LIBRARY=/usr/lib64/liblapack.a \
#  -D HDF5_LIBRARIES=/home/yf30/alps-2.2.b3-r7462-src-with-boost/hdf5_installed/lib/libhdf5.a \
#  -D CMAKE_INSTALL_PREFIX=/home/yf30/alps_static_nano/ \
#  -D LAPACK_64_BIT=ON \
#  -D HDF5_INCLUDE_DIR=/home/yf30/alps-2.2.b3-r7462-src-with-boost/hdf5_installed/include/ \
#  -D Boost_ROOT_DIR=/home/yf30/alps-2.2.b3-r7462-src-with-boost/boost/ \
#  -D BUILD_SHARED_LIBS=OFF \
#  -D PYTHON_LIBRARY=/usr/lib64/libpython2.6.so \
#  -D CMAKE_C_COMPILER=/usr/bin/gcc \
#  -D CMAKE_CXX_COMPILER=/usr/bin/g++ \
#  -D ALPS_ENABLE_MPI=ON \
#  -D MPI_CXX_COMPILER=/usr/lib64/openmpi/bin/mpicxx \
#  -D MPI_C_INCLUDE_PATH=/usr/include/openmpi-x86_64/ \
#  -D MPI_C_LIBRARIES="/usr/lib64/openmpi/lib/libmpi.so;/usr/lib64/openmpi/lib/libopen-pal.so;/usr/lib64/openmpi/lib/libopen-rte.so" \
#  -D MPI_CXX_INCLUDE_PATH=/usr/include/openmpi-x86_64/ \
#  -D MPI_CXX_LIBRARIES="/usr/lib64/openmpi/lib/libmpi_cxx.so;/usr/lib64/openmpi/lib/libopen-pal.so;/usr/lib64/openmpi/lib/libopen-rte.so" \
#  /home/yf30/alps-2.2.b3-r7462-src-with-boost/alps/
#
##  End of working version

elif [ "$1" == "shared" ]; then

  echo -e "cmake: build a shared (dynamical) ALPS library..."
  cmake \
  -D ALPS_BUILD_TESTS=OFF \
  -D ALPS_BUILD_EXAMPLES=OFF \
  -D ALPS_BUILD_APPLICATIONS=ON \
  -D ALPS_BUILD_PYTHON=ON \
  -D CMAKE_INSTALL_PREFIX=/home/yf30/alps/ \
  -D LAPACK_64_BIT=ON \
  -D Boost_ROOT_DIR=/home/yf30/alps-2.2.b3-r7462-src-with-boost/boost/ \
  -D CMAKE_C_COMPILER=/usr/bin/gcc \
  -D CMAKE_CXX_COMPILER=/usr/bin/g++ \
  -D ALPS_ENABLE_MPI=ON \
  /home/yf30/alps-2.2.b3-r7462-src-with-boost/alps/

else

  echo -e "To build the ALPS library, this script takes one argument (\"static\" or \"shared\")."

fi

############################ working area ############################

#  -D BLAS_LIBRARY=/usr/lib64/libblas.so \
#  -D LAPACK_LIBRARY=/usr/lib64/liblapack.so \
#  -D MPI_INCLUDE_PATH=/usr/include/openmpi-x86_64/ \
#  -D MPI_CXX_COMPILER=/usr/lib64/openmpi/bin/mpicxx \
#  -D MPI_C_LIBRARIES="/usr/lib64/openmpi/lib/libmpi.so;/usr/lib64/openmpi/lib/libopen-pal.so;/usr/lib64/openmpi/lib/libopen-rte.so" \
#  -D MPI_CXX_LIBRARIES="/usr/lib64/openmpi/lib/libmpi_cxx.so;/usr/lib64/openmpi/lib/libopen-pal.so;/usr/lib64/openmpi/lib/libopen-rte.so" \


#  cmake \
#  -D ALPS_BUILD_TESTS=OFF \
#  -D ALPS_BUILD_EXAMPLES=OFF \
#  -D ALPS_BUILD_APPLICATIONS=ON \
#  -D ALPS_BUILD_PYTHON=ON \
#  -D BLAS_LIBRARY=/usr/lib64/libblas.so \
#  -D LAPACK_LIBRARY=/usr/lib64/liblapack.so \
#  -D CMAKE_INSTALL_PREFIX=/home/yf30/alps/ \
#  -D LAPACK_64_BIT=ON \
#  -D Boost_ROOT_DIR=/home/yf30/alps-2.2.b3-r7462-src-with-boost/boost/ \
#  -D CMAKE_C_COMPILER=/usr/bin/gcc \
#  -D CMAKE_CXX_COMPILER=/usr/bin/g++ \
#  -D ALPS_ENABLE_MPI=ON \
#  -D MPI_CXX_COMPILER=/usr/lib64/openmpi/bin/mpicxx \
#  -D MPI_C_INCLUDE_PATH=/usr/include/openmpi-x86_64/ \
#  -D MPI_C_LIBRARIES="/usr/lib64/openmpi/lib/libmpi.so;/usr/lib64/openmpi/lib/libopen-pal.so;/usr/lib64/openmpi/lib/libopen-rte.so" \
#  -D MPI_CXX_INCLUDE_PATH=/usr/include/openmpi-x86_64/ \
#  -D MPI_CXX_LIBRARIES="/usr/lib64/openmpi/lib/libmpi_cxx.so;/usr/lib64/openmpi/lib/libopen-pal.so;/usr/lib64/openmpi/lib/libopen-rte.so" \
#  /home/yf30/alps-2.2.b3-r7462-src-with-boost/alps/

# 
#  -D MPI_LIBRARIES="/usr/lib64/openmpi/lib/libmpi.so;/usr/lib64/openmpi/lib/libmpi_cxx.so;/usr/lib64/openmpi/lib/libopen-pal.so;/usr/lib64/openmpi/lib/libopen-rte.so" \
#  -D HDF5_LIBRARIES=/usr/lib64/libhdf5.so.6 \
#  -D HDF5_INCLUDE_DIR=/home/yf30/alps-2.2.b3-r7462-src-with-boost/hdf5_installed/include/ \
#  -D PYTHON_LIBRARY=/usr/lib64/libpython2.6.so \
#-D CMAKE_EXE_LINKER_FLAGS="${CMAKE_EXE_LINKER_FLAGS} -lgfortran" \
#-D CMAKE_EXE_LINKER_FLAGS=-lgfortran \
