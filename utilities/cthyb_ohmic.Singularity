BootStrap: yum
OSVersion: 7
MirrorURL: http://mirror.centos.org/centos-%{OSVERSION}/%{OSVERSION}/os/$basearch/
Include: yum

%setup
# get alps
wget http://alps.comp-phys.org/static/software/releases/alps-2.2.b3-r7462-src-with-boost.tar.gz
tar -xzvf alps-2.2.b3-r7462-src-with-boost.tar.gz 1> /dev/null
mv alps-2.2.b3-r7462-src-with-boost/ ${SINGULARITY_ROOTFS}/opt/
rm alps-2.2.b3-r7462-src-with-boost.tar.gz

# get cthyb-ohmic
git clone -b color_swap https://github.com/leofang/cthyb-ohmic.git
mv cthyb-ohmic/ ${SINGULARITY_ROOTFS}/opt/

%environment 
PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" 
export PATH

%runscript 
exec /usr/local/bin/cthyb_ohmic "$@"

%post
# install build tools
yum -y install epel-release # for hdf5
yum -y install gcc gcc-c++ gcc-gfortran make cmake blas blas-devel lapack lapack-devel hdf5 hdf5-devel

# install alps
mkdir /opt/alps_build
cd /opt/alps_build
cmake \
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
make -j 2
make install

# install cthyb-ohmic
cd ../cthyb-ohmic/
cmake ./
make -j 2
mv cthyb_ohmic /usr/local/bin/

# clean up
cd ..
rm -rf cthyb-ohmic/ alps-2.2.b3-r7462-src-with-boost/ alps_build/
yum clean all

%labels
Maintainer  Leo Fang <leofang@phy.duke.edu>
Description A Singularity container for cthyb-ohimc
License     WTFPL v2
