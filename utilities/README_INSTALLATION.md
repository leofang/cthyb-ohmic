This instruction is for intalling static version of ALPS and cthyb_ohmic on 
a Physics machine, and is tested with the following environment:

g++ v4.4.7; HDF5 1.8.12; ALPS v.2.2.b3 (with Boost 1.54); 

test machine: lqcd28 & nano06 (other machines may also work if Jimmy does not 
mess up the installed packages...)

Last updated: Dec. 6, 2015


Install the ALPS library:

1. Download ALPS and extract the files to a source folder, say ~/alps_boost_src/
   There will be two subfolders (alps and boost) in this directory.

2. If you are building the library on nano machines (nano06~09), skip the HDF5
   part and go to step 8 because a static version of HDF5 is installed on these
   machines.

3. Download HDF5 and extract the files to a source folder, say ~/hdf5_src/.

4. Create a build folder for HDF5 (eg. ~/hdf5_build/), and go to the folder.

5. Point the key "CMAKE_INSTALL_PREFIX" in the script "cmake_input_hdf5.sh" to 
   ~/alps_boost_src/hdf5/ so that the HDF5 library will be installed there.
   Also modify the last line of the script (path to the HDF5 source).

6. Execute the script "cmake_input_hdf5" in the HDF5 build folder

7. "make" and then "make install"

8. Create a build folder for ALPS (eg. ~/alps_build/), and go to the folder.

9. Modify the script "cmake_input.sh" accordingly (let's install it to ~/alps/) 
   and execute it in the ALPS build forder.

10. Repeat step 7.

11. Done. Now one can check that in ~/alps/lib/ there are two major static 
    libraries: libalps.a and libboost.a.


Install the cthyb_ohmic application:

1. Download its source code, extract it to ~/cthyb/, and go to this folder

2. Add a line to CMakeLists.txt before the "find_package" line so that cmake 
   can find the ALPS library: set(ALPS_ROOT_DIR "~/alps")

3. Type cmake ./ and then make 

4. Done, an executable "cthyb_ohmic" is placed in the same folder.


Note:

1. The Python support for the ALPS library (pyalps) will not be installed
   because we build static libraries to which Python cannot dynamically link.
   To access this utility, one has to build a shared version of ALPS and 
   use alpspython therein. This is supported by the script cmake_input.sh
   with argument "shared."

2. Strictly speaking, any ALPS application like cthyb_ohmic built here is 
   NOT fully static, as one can check using ldd that it is still dynamically
   linked to system libraries (such as libstdc++.so). This is the best we
   can do because most of these libraries do not have a static version.

3. To install on workstations over which you have some control (such as 
   running package managers like pkcon), use the script "cmake_input_grads-bd.sh"
   instead of "cmake_input.sh".


Known issues:

1. The source code of alps/hybridization has a bug, see https://github.com/ALPSCore/alps-cthyb/issues/4.

2. The mpi.py script may not work (on some machines it complains "cannot find 
   boost.mpi"). Change "boost.mpi" to "mpi" would fix it. 


Useful links:

https://alps.comp-phys.org/mediawiki/index.php/Building_ALPS_from_source

https://alps.comp-phys.org/mediawiki/index.php/ALPS2_cmake_options

https://alps.comp-phys.org/trac/wiki/AlpsStatic

https://lists.phys.ethz.ch/pipermail/comp-phys-alps-users/2014/002314.html
