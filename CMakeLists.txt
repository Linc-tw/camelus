

    ######################################################
    ##  CMakeLists.txt for Camelus v2.0			##
    ##  Version 2018.03.05				##
    ##							##
    ##  Copyright (C) 2018 - Chieh-An Lin		##
    ##  GNU GPLv3 - https://www.gnu.org/licenses/	##
    ######################################################


cmake_minimum_required(VERSION 2.6)
project(camelus)

set(CMAKE_BUILD_TYPE Debug)

## Use pkg-config to find modules
find_package(PkgConfig REQUIRED)
pkg_check_modules(REQ_PKG REQUIRED gsl fftw3)
pkg_check_modules(CFITSIO cfitsio)
pkg_check_modules(CHEALPIX chealpix)
pkg_check_modules(HEALPIX_CXX healpix_cxx)
pkg_check_modules(MPI mpi)

## Link the required modules
include_directories(${REQ_PKG_INCLUDE_DIRS})
link_directories(${REQ_PKG_LIBRARY_DIRS})
include_directories($ENV{NICAEA}/include)
link_directories($ENV{NICAEA}/lib)
set(L_FLAGS "${REQ_PKG_LIBRARIES} ${REQ_PKG_INCLUDE} -lnicaea")

## Optional module cfitsio, link if found
if (CFITSIO_FOUND)
  add_definitions(-D __CAMELUS_USE_FITS__)
  include_directories(${CFITSIO_INCLUDE_DIRS})
  link_directories(${CFITSIO_LIBRARY_DIRS})
endif()

## Optional module chealpix, link if found
if (CHEALPIX_FOUND)
  add_definitions(-D __CAMELUS_USE_HEALPIX__)
  include_directories(${CHEALPIX_INCLUDE_DIRS})
  link_directories(${CHEALPIX_LIBRARY_DIRS})
endif()

## Optional module healpix_cxx, link if found
if (HEALPIX_CXX_FOUND)
  add_definitions(-D __CAMELUS_USE_HEALPIX_CXX__)
  include_directories(${HEALPIX_CXX_INCLUDE_DIRS})
  link_directories(${HEALPIX_CXX_LIBRARY_DIRS})
endif()

## Optional module mpi, link if found
if (MPI_FOUND)
  add_definitions(-D __CAMELUS_USE_MPI__)
  include_directories(${MPI_INCLUDE_DIRS})
  link_directories(${MPI_LIBRARY_DIRS})
endif()

## Define cflags
set(DEF_C_FLAGS "-Wall -Wpedantic -Wno-unused-result -Wno-maybe-uninitialized -O3 -ffast-math")
set(ADD_C_FLAGS "-Wno-unused-variable")
set(CMAKE_C_FLAGS "${DEF_C_FLAGS} ${ADD_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${DEF_C_FLAGS} ${ADD_C_FLAGS}")

## Build executable
#file(GLOB_RECURSE SRC_FILES "${PROJECT_SOURCE_DIR}/source/*.c" "${PROJECT_SOURCE_DIR}/source/*.cpp")
file(GLOB_RECURSE SRC_FILES "${PROJECT_SOURCE_DIR}/source/*.c")
if (CHEALPIX_FOUND)
  file(GLOB_RECURSE HEALPIX_FILES "${PROJECT_SOURCE_DIR}/source/wrapper.cpp")
endif()
add_executable(camelus ${SRC_FILES} ${HEALPIX_FILES})
target_link_libraries(camelus
  ${REQ_PKG_LIBRARIES}     ${REQ_PKG_INCLUDE}
  -lnicaea
  ${CFITSIO_LIBRARIES}     ${CFITSIO_INCLUDE}
  ${CHEALPIX_LIBRARIES}    ${CHEALPIX_INCLUDE}
  ${HEALPIX_CXX_LIBRARIES} ${HEALPIX_CXX_INCLUDE}
  ${MPI_LIBRARIES}         ${MPI_INCLUDE}
)

