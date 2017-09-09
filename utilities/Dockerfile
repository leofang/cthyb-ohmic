FROM centos:latest

# get alps
#RUN yum -y install wget git
ADD http://alps.comp-phys.org/static/software/releases/alps-2.2.b3-r7462-src-with-boost.tar.gz /opt/
#RUN tar -xzvf alps-2.2.b3-r7462-src-with-boost.tar.gz -C /opt/ 1> /dev/null

# get cthyb-ohmic; note the relative path
ADD . /opt/cthyb-ohmic/

# check everything is in place
RUN ls -l /opt/
RUN ls -l /opt/cthyb-ohmic/
RUN ls -l /opt/alps-2.2.b3-r7462-src-with-boost/

# install build tools
RUN yum -y install epel-release # for hdf5
RUN yum -y install gcc gcc-c++ gcc-gfortran make cmake blas blas-devel lapack lapack-devel hdf5 hdf5-devel

# install alps
RUN mkdir /opt/alps_build
RUN cd /opt/alps_build
RUN cmake \
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
      /opt/alps-2.2.b3-r7462-src-with-boost/alps/;
RUN make -j 2
RUN make install

# install cthyb-ohmic
RUN cd /opt/cthyb-ohmic/
RUN cmake ./
RUN make -j 2
RUN mv cthyb_ohmic /usr/local/bin/

# clean up
RUN yum clean all
RUN cd /opt/
RUN rm -rf cthyb-ohmic/ alps-2.2.b3-r7462-src-with-boost/ alps_build/ alps-2.2.b3-r7462-src-with-boost.tar.gz

# catch current time
ARG build_time=$(date +"%b. %d, %T %Z, %Y")

LABEL Maintainer  = "Leo Fang <leofang@phy.duke.edu>" \
      Description = "A Singularity container for cthyb-ohmic" \
      License     = "WTFPL v2" \
      Build_date  = $build_time
