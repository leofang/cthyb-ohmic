#  Copyright Matthias Troyer 2009 - 2013.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(alpsize NONE)

message("cmake starts...")

#find_package(blas REQUIRED)
#message(STATUS "Found BLAS")
#
#find_package(lapack REQUIRED)
#message(STATUS "Found Lapack")

# find ALPS Library
find_package(ALPS REQUIRED PATHS ${ALPS_ROOT_DIR} $ENV{ALPS_HOME} NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "Found ALPS: ${ALPS_ROOT_DIR} (revision: ${ALPS_VERSION})")
#message("ALPS_USE_FILE is ${ALPS_USE_FILE}")
include(${ALPS_USE_FILE})

# BLAS/LAPACK
#message("BLAS_FOUND: ${ALPS_BLAS_FOUND}")
#message("LAPACK_FOUND: ${ALPS_LAPACK_FOUND}")
#message("LAPACK_DEFINITIONS: ${ALPS_LAPACK_DEFINITIONS}")
#message("LAPACK_LINKER_FLAGS: ${ALPS_LAPACK_LINKER_FLAGS}")
#message("LAPACK_LIBRARIES: ${ALPS_LAPACK_LIBRARIES}")
#message("LAPACK_LIBRARY: ${ALPS_LAPACK_LIBRARY}")
#message("BLAS_LIBRARIES: ${ALPS_BLAS_LIBRARIES}")
#message("BLAS_LIBRARY: ${ALPS_BLAS_LIBRARY}")
#message("MKL_INCLUDE_DIR: ${ALPS_MKL_INCLUDE_DIR}")

# enable C++ compilers
enable_language(C CXX)

# rule for generating 'hello world' program
#add_executable(hello hello.C)
#target_link_libraries(hello ${ALPS_LIBRARIES})


if(LAPACK_FOUND)

#  set (DMFT_QMC_SOURCE
#       externalsolver.C 
#       fouriertransform.C 
#       hilberttransformer.C 
#       hirschfyesim.C 
#       alps_solver.C 
#       bandstructure.C 
#       auxiliaryfunctions.C)
#  add_library(dmft_qmc_impl STATIC ${DMFT_QMC_SOURCE})

#  set (DMFT_QMC_INTERACTION_EXPANSION_SOURCE 
#       interaction_expansion/auxiliary.cpp
#       interaction_expansion/observables.cpp 
#       interaction_expansion/fastupdate.cpp
#       interaction_expansion/selfenergy.cpp 
#       interaction_expansion/solver.cpp 
#       interaction_expansion/io.cpp 
#       interaction_expansion/splines.cpp 
#       interaction_expansion/interaction_expansion.cpp 
#       interaction_expansion/measurements.cpp 
#      interaction_expansion/model.cpp)
#
#  set (DMFT_QMC_INTERACTION_EXPANSION2_SOURCE
#       fouriertransform.C
#       interaction_expansion2/auxiliary.cpp
#       interaction_expansion2/observables.cpp 
#       interaction_expansion2/fastupdate.cpp
#       interaction_expansion2/selfenergy.cpp 
#       interaction_expansion2/solver.cpp 
#       interaction_expansion2/io.cpp 
#       interaction_expansion2/splines.cpp 
#       interaction_expansion2/interaction_expansion.cpp 
#       interaction_expansion2/measurements.cpp 
#       interaction_expansion2/model.cpp)

  set (DMFT_QMC_HYBRIDIZATION_SOURCE
       hybsim.cpp
       hyblocal.cpp
       hybint.cpp
       hybfun.cpp
       hybretintfun.cpp
       hybmatrix.cpp
       hybmatrix_ft.cpp
       hybconfig.cpp
       hybupdates.cpp
       hybevaluate.cpp
       hybmeasurements.cpp)

#  set (DMFT_QMC_SOLVERS "")
 
  #Interaction Expansion
#  set (DMFT_QMC_SOLVERS ${DMFT_QMC_SOLVERS} interaction_expansion_impl)
#  add_library (interaction_expansion_impl STATIC ${DMFT_QMC_INTERACTION_EXPANSION_SOURCE})

#  target_link_libraries(dmft_qmc_impl alps ${LAPACK_LIBRARY} ${BLAS_LIBRARY} ${ALPS_LIBRARIES})
  add_definitions(${LAPACK_DEFINITIONS} )
#  include_directories(${PROJECT_SOURCE_DIR}/applications/dmft/qmc)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libstdc++")

#  add_executable(dmft main.C selfconsistency.C)
#  add_executable(hirschfye solver_main.C hirschfyesim.C)

  add_executable(hybridization hybmain.cpp ${DMFT_QMC_HYBRIDIZATION_SOURCE})
#  add_executable(interaction interaction_expansion2/main.cpp ${DMFT_QMC_INTERACTION_EXPANSION2_SOURCE})

#  target_link_libraries(dmft dmft_qmc_impl  ${DMFT_QMC_SOLVERS} ${ALPS_LIBRARIES})
#  target_link_libraries(hirschfye dmft_qmc_impl ${ALPS_LIBRARIES})

  target_link_libraries(hybridization ${LAPACK_LIBRARY} ${BLAS_LIBRARY} ${ALPS_LIBRARIES})
#  target_link_libraries(hybridization dmft_qmc_impl ${ALPS_LIBRARIES})
#  target_link_libraries(interaction dmft_qmc_impl)

#  install(TARGETS dmft RUNTIME DESTINATION bin COMPONENT applications)
#  install(TARGETS hirschfye  RUNTIME DESTINATION bin COMPONENT applications)
#  install(TARGETS hybridization RUNTIME DESTINATION bin COMPONENT applications)
#  install(TARGETS interaction RUNTIME DESTINATION bin COMPONENT applications)
#  install(FILES hybridization/Documentation/hybdoc.pdf DESTINATION doc)

#  message("LAPACK_DEFINITIONS=${LAPACK_DEFINITIONS}")
#  message("${CMAKE_EXE_LINKER_FLAGS}") 
#  message("${LAPACK_LINKER_FLAGS}")
#  message("${LAPACK_LIBRARY}") 
#  message("${BLAS_LIBRARY}") 
#  message("${ALPS_LIBRARIES}")

else(LAPACK_FOUND)
  message(STATUS "dmft will not be built since the lapack library has not been found")
endif(LAPACK_FOUND)
