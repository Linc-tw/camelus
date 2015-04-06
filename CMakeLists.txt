

    ##############################################
    ##  CMakeLists.txt for camelus v1.2		##
    ##  Chieh-An Lin				##
    ##  Version 2015.04.05			##
    ##############################################

  
cmake_minimum_required(VERSION 2.6)
include(FindPkgConfig)

project(camelus)

pkg_check_modules(PKGS REQUIRED gsl fftw3)
include_directories(${PKGS_INCLUDE_DIRS})
link_directories(${PKGS_LIBRARY_DIRS})

include_directories($ENV{NICAEA}/include)
link_directories($ENV{NICAEA}/lib)

set(CMAKE_C_FLAGS "-O3 -ffast-math")
set(CMAKE_CXX_FLAGS "-O3 -ffast-math")

FILE(GLOB_RECURSE src_files "${PROJECT_SOURCE_DIR}/source/*.c")
add_executable(camelus ${src_files})
target_link_libraries(camelus -lnicaea ${PKGS_LIBRARIES} ${PKGS_INCLUDE})