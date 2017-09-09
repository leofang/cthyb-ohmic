#!/bin/bash
# Leo Fang <leofang@phy.duke.edu>
# Create a Singularity image for the cthyb-ohmic application
# Last update: Sep 8, 2017


image_filename="cthyb_ohmic" # creating an image named ${image_filename}.img
temp_dir="/opt/"             # temporary directory in the container 

# prepare a CentOS7 image and install build tools
# TODO: mark specific versions for better reproducibility?
singularity create ${image_filename}.img
singularity import ${image_filename}.img docker://centos:latest
sudo singularity exec --writable ${image_filename}.img \
    yum -y install epel-release # for hdf5
sudo singularity exec --writable ${image_filename}.img \
    yum -y install gcc gcc-c++ gcc-gfortran make cmake blas blas-devel lapack lapack-devel hdf5 hdf5-devel

# get alps
wget http://alps.comp-phys.org/static/software/releases/alps-2.2.b3-r7462-src-with-boost.tar.gz
tar -xzvf alps-2.2.b3-r7462-src-with-boost.tar.gz

# write to the host's space, so no need to sudo
mkdir ${image_filename}_build
singularity exec -B ./:${temp_dir} --pwd ${temp_dir}/${image_filename}_build --writable ${image_filename}.img bash -c \
     "cmake \
        -D ALPS_BUILD_TESTS=OFF \
        -D ALPS_BUILD_EXAMPLES=OFF \
        -D ALPS_BUILD_APPLICATIONS=OFF \
        -D ALPS_BUILD_PYTHON=OFF \
        -D ALPS_ENABLE_MPI=OFF \
        -D LAPACK_64_BIT=ON \
        -D CMAKE_C_COMPILER=/usr/bin/gcc \
        -D CMAKE_CXX_COMPILER=/usr/bin/g++ \
        -D CMAKE_INSTALL_PREFIX=/usr/local/ \
        -D Boost_ROOT_DIR=${temp_dir}/alps-2.2.b3-r7462-src-with-boost/boost/ \
        ${temp_dir}/alps-2.2.b3-r7462-src-with-boost/alps/;
      make -j 2 "

# now sudo to install files to /usr/local
sudo singularity exec -B ./:${temp_dir} --pwd ${temp_dir}/${image_filename}_build --writable ${image_filename}.img make install

# clone my code
git clone -b color_swap https://github.com/leofang/cthyb-ohmic.git

# install cthyb with sudo
sudo singularity exec -B ./:${temp_dir} --pwd ${temp_dir}/cthyb-ohmic --writable ${image_filename}.img \
    bash -c "cmake -D CMAKE_INSTALL_PREFIX=/usr/local/bin/ ./; make -j 2; mv cthyb_ohmic /usr/local/bin"

# clean up
sudo rm -rf alps-2.2.b3-r7462-src-with-boost.tar.gz alps-2.2.b3-r7462-src-with-boost/ cthyb-ohmic/ ${image_filename}_build/

# modify runscript 
echo -e '#!/bin/bash\nexec /usr/local/bin/cthyb_ohmic "$@"' > runscript
chmod +x runscript
sudo singularity exec -B ./:${temp_dir} --writable ${image_filename}.img \
    bash -c "chown root:root ${temp_dir}/runscript; mv ${temp_dir}/runscript /.singularity.d/;"

# modify metadata
echo "{" > labels.json
echo "   \"Maintainer\"  : \"Leo Fang <leofang@phy.duke.edu>\"" >> labels.json
echo "   \"Build_date\"  : \"$(date +"%b. %d, %T %Z, %Y")\"" >> labels.json 
echo "   \"Description\" : \"A Singularity container for cthyb-ohimc\"" >> labels.json
echo "   \"License\"     : \"WTFPL v2\"" >> labels.json
echo "}" >> labels.json
chmod +x labels.json
sudo singularity exec -B ./:${temp_dir} --writable ${image_filename}.img \
    bash -c "chown root:root ${temp_dir}/labels.json; mv ${temp_dir}/labels.json /.singularity.d/;"
