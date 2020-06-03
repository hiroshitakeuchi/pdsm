cmake_minimum_required(VERSION 3.16)
project(pdsm)

set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

include_directories(
        /usr/local/include/eigen3/
#        /usr/local/include/boost/
#        /usr/local/include/boost/functional/
#        ./include
)

add_executable(pdsm main.cpp include/pdsm.h include/rips.h include/persistent_homology.h include/common_definitions_pdsm.h)

add_executable(pdsm_torus main_torus.cpp include/pdsm.h include/rips_torus.h include/persistent_homology.h include/common_definitions_pdsm.h)
