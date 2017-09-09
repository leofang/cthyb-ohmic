FROM centos:latest

# get cthyb-ohmic
ADD . /opt/cthyb-ohmic/

# get alps
RUN yum -y install wget \
    && wget -q http://alps.comp-phys.org/static/software/releases/alps-2.2.b3-r7462-src-with-boost.tar.gz \
    && tar -xzf alps-2.2.b3-r7462-src-with-boost.tar.gz -C /opt/

# install build tools
RUN yum -y install epel-release \
    && yum -y install gcc gcc-c++ gcc-gfortran make cmake blas blas-devel lapack lapack-devel hdf5 hdf5-devel

# install alps
RUN mkdir /opt/alps_build && cd /opt/alps_build && cmake \
      -D ALPS_BUILD_TESTS=OFF \
      -D ALPS_BUILD_EXAMPLES=OFF \
      -D ALPS_BUILD_APPLICATIONS=OFF \
      -D ALPS_BUILD_PYTHON=OFF \
      -D ALPS_ENABLE_MPI=OFF \
      -D LAPACK_64_BIT=ON \
      -D CMAKE_C_COMPILER=/usr/bin/gcc \
      -D CMAKE_CXX_COMPILER=/usr/bin/g++ \
      -D CMAKE_INSTALL_PREFIX=/usr/local/ \
      -D Boost_ROOT_DIR=/opt/alps-2.2.b3-r7462-src-with-boost/boost/ \
      /opt/alps-2.2.b3-r7462-src-with-boost/alps/ \
    && make && make install

# install cthyb-ohmic
RUN cd /opt/cthyb-ohmic/ && cmake ./ && make && mv cthyb_ohmic /usr/local/bin/

# clean up
RUN yum clean all && rm -rf /opt/

# make a label
# TODO: find a way to label the build-complete time
LABEL Maintainer="Leo Fang <leofang@phy.duke.edu>" Description="A Singularity container for cthyb-ohmic" License="WTFPL v2"
