cmake_minimum_required (VERSION 3.0)


project(allele_age_simulator)

if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_VERBOSE_MAKEFILE ON)

if(LINUX) 
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ~/.linuxbrew/share/cmake/Modules)
    list(APPEND CMAKE_PREFIX_PATH "~/.linuxbrew")
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -Wall")

find_package(Armadillo REQUIRED)
find_package(GSL REQUIRED)
find_package(OpenMP REQUIRED)


include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${GSL_INCLUDE_DIRS})

add_executable(allele_age_simulator allele_age_simulator.cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
target_link_libraries(allele_age_simulator ${ARMADILLO_LIBRARIES})
target_link_libraries(allele_age_simulator ${GSL_LIBRARIES})

