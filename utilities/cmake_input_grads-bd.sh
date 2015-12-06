#!/bin/bash
# tweaked for grads-bd

echo -e
echo -e "This script is designed for installing ALPS on my workstation (grads-bd) only!"
echo -e

if [ "$1" == "static" ]; then

  # My gcc version is 5.1.1 which is too new for the Boost source code (v1.54) that comes
  # with the ALPS library. Therefore, ALPS has to be compiled against the newer version 
  # (v1.57) installed using pkcon on my computer.
  #
  # The extra linker flags are necessary after a lot of experiments (they are here mainly
  # because of the static version of HDF5 installed by pkcon), so do not remove them!

  echo -e "cmake: build a static ALPS library..."
  echo -e "After running cmake, modify libhdf5.so and libhdf5_hl.so in CMakeCache.txt to their"
  echo -e "static counterparts."
  echo -e
  cmake \
  -D CMAKE_INSTALL_PREFIX=/home/yf30/alps_workstation_static/ \
  -D CMAKE_EXE_LINKER_FLAGS="${CMAKE_EXE_LINKER_FLAGS} -lpthread -ldl -lc -lgfortran" \
  -D CMAKE_FIND_LIBRARY_SUFFIXES=".a" \
  -D BLAS_LIBRARY=/usr/lib64/libblas.a \
  -D LAPACK_LIBRARY=/usr/lib64/liblapack.a \
  -D LAPACK_64_BIT=ON \
  -D BUILD_SHARED_LIBS=OFF \
  -D Boost_USE_STATIC_LIBS=ON \
  -D HDF5_LIBRARIES="/usr/lib64/libhdf5.a;/usr/lib64/libhdf5_hl.a" \
  -D ALPS_BUILD_TESTS=OFF \
  -D ALPS_BUILD_EXAMPLES=OFF \
  -D ALPS_BUILD_APPLICATIONS=ON \
  -D ALPS_BUILD_PYTHON=OFF \
  -D CMAKE_C_COMPILER=/usr/bin/gcc \
  -D CMAKE_CXX_COMPILER=/usr/bin/g++ \
  -D ALPS_ENABLE_MPI=ON \
  /home/yf30/alps-2.2.b3-r7462-src-with-boost/alps/

elif [ "$1" == "shared" ]; then

  echo -e "cmake: build a shared (dynamical) ALPS library..."
  cmake \
  -D ALPS_BUILD_TESTS=OFF \
  -D ALPS_BUILD_EXAMPLES=OFF \
  -D ALPS_BUILD_APPLICATIONS=ON \
  -D ALPS_BUILD_PYTHON=ON \
  -D CMAKE_INSTALL_PREFIX=/home/yf30/alps_workstation/ \
  -D LAPACK_64_BIT=ON \
  -D CMAKE_C_COMPILER=/usr/bin/gcc \
  -D CMAKE_CXX_COMPILER=/usr/bin/g++ \
  -D ALPS_ENABLE_MPI=ON \
  /home/yf30/alps-2.2.b3-r7462-src-with-boost/alps/

else

  echo -e "To build the ALPS library, this script takes one argument (\"static\" or \"shared\")."

fi


######################## working area ########################


#  -D HDF5_USE_STATIC_LIBRARIES=ON \
#  -D CMAKE_EXE_LINKER_FLAGS="${CMAKE_EXE_LINKER_FLAGS} -lpthread -ldl -lc" \

#  -D CMAKE_EXE_LINKER_FLAGS="-lboost_mpi -lboost_serialization" \
#  -D PYTHON_LIBRARY="/usr/lib64/libpython2.7.so;/usr/lib64/openmpi/lib/libboost_graph_parallel.so;/usr/lib64/openmpi/lib/libboost_mpi.so;/usr/lib64/openmpi/lib/libboost_mpi_python.so;/usr/lib64/openmpi/lib/mpi.so;/usr/lib64/openmpi/lib/libmpi.so;/usr/lib64/openmpi/lib/libmpi_cxx.so;/usr/lib64/openmpi/lib/libopen-pal.so;/usr/lib64/openmpi/lib/libopen-rte.so" \
#  -D Boost_MPI_LIBRARY="/usr/lib64/openmpi/lib/libboost_graph_parallel.so;/usr/lib64/openmpi/lib/libboost_mpi.so;/usr/lib64/openmpi/lib/libboost_mpi_python.so;/usr/lib64/openmpi/lib/mpi.so" \
#  -D CMAKE_CXX_FLAGS=-std=gnu++98 \
#  -D BLAS_LIBRARY=/usr/lib64/libblas.a \
#  -D LAPACK_LIBRARY=/usr/lib64/liblapack.a \
#  -D HDF5_LIBRARIES=/usr/lib64/libhdf5.a \
#  -D CMAKE_INSTALL_PREFIX=/home/yf30/alps_static/ \
#  -D LAPACK_64_BIT=ON \
#  -D HDF5_INCLUDE_DIR=/usr/include/ \
#  -D Boost_ROOT_DIR=/home/yf30/alps-2.2.b3-r7462-src-with-boost/boost/ \
#  -D BUILD_SHARED_LIBS=OFF \
#  -D ALPS_BUILD_APPLICATIONS=OFF \
#  -D PYTHON_LIBRARY=/usr/lib64/libpython2.7.so \
#  -D CMAKE_C_COMPILER=/usr/bin/gcc \
#  -D CMAKE_CXX_COMPILER=/usr/bin/g++ \
#  -D ALPS_ENABLE_MPI=ON \
#  -D MPI_CXX_COMPILER=/usr/lib64/openmpi/bin/mpicxx \
#  -D MPI_C_INCLUDE_PATH=/usr/include/openmpi-x86_64/ \
#  -D MPI_C_LIBRARIES=/usr/lib64/openmpi/lib/libmpi.so \
#  -D MPI_CXX_INCLUDE_PATH=/usr/include/openmpi-x86_64/ \
#  -D MPI_CXX_LIBRARIES=/usr/lib64/openmpi/lib/libmpi_cxx.so \
#  /home/yf30/alps-2.2.b3-r7462-src-with-boost/alps/


#  -D Boost_MPI_LIBRARY="/usr/lib64/openmpi/lib/libboost_mpi.so.1.57.0;/usr/lib64/openmpi/lib/libboost_mpi_python.so.1.57.0;/usr/lib64/openmpi/lib/mpi.so" \
#  -D HDF5_LIBRARIES=/usr/lib64/libhdf5.so \
#  -D HDF5_INCLUDE_DIR=/usr/include \
#  -D PYTHON_LIBRARY=/usr/lib64/libpython2.7.so \
#  -D BLAS_LIBRARY=/usr/lib64/libblas.so \
#  -D LAPACK_LIBRARY=/usr/lib64/liblapack.so \
#  -D MPIEXEC=/usr/lib64/openmpi/bin/mpiexec \
#  -D MPI_CXX_INCLUDE_PATH=/usr/include/openmpi-x86_64/ \
#  -D MPI_CXX_COMPILER=/usr/lib64/openmpi/bin/mpicxx \
#  -D MPI_C_INCLUDE_PATH=/usr/include/openmpi-x86_64/ \
#  -D MPI_C_LIBRARIES="/usr/lib64/openmpi/lib/libmpi.so;/usr/lib64/openmpi/lib/libopen-pal.so;/usr/lib64/openmpi/lib/libopen-rte.so;" \
#  -D MPI_CXX_LIBRARIES="/usr/lib64/openmpi/lib/libmpi.so;/usr/lib64/openmpi/lib/libmpi_cxx.so;/usr/lib64/openmpi/lib/libopen-pal.so;/usr/lib64/openmpi/lib/libopen-rte.so;" \




#  -D CMAKE_BUILD_TYPE=Debug \

#  -D Boost_LIBRARY_DIR=/usr/lib64/ \
#  -D MPI_C_INCLUDE_PATH=/usr/include/openmpi-x86_64/ \
#  -D MPI_C_LIBRARIES="/usr/lib64/openmpi/lib/libmpi.so;/usr/lib64/openmpi/lib/libopen-pal.so;/usr/lib64/openmpi/lib/libopen-rte.so;/usr/lib64/openmpi/lib/libboost_mpi.so.1.57.0" \
#  -D MPI_CXX_INCLUDE_PATH=/usr/include/openmpi-x86_64/ \
#  -D MPI_CXX_LIBRARIES="/usr/lib64/openmpi/lib/libmpi.so;/usr/lib64/openmpi/lib/libmpi_cxx.so;/usr/lib64/openmpi/lib/libopen-pal.so;/usr/lib64/openmpi/lib/libopen-rte.so;/usr/lib64/openmpi/lib/libboost_mpi.so.1.57.0" \


#-D Boost_ROOT_DIR=/home/yf30/alps-2.2.b3-r7462-src-with-boost/boost/ \
#-D CMAKE_EXE_LINKER_FLAGS="${CMAKE_EXE_LINKER_FLAGS} -lgfortran" \
#-D CMAKE_EXE_LINKER_FLAGS=-lgfortran \
