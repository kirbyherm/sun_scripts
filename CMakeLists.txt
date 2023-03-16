
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)

#set(CMAKE_VERBOSE_MAKEFILE ON) 

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules/")
message(STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH}) 
project(gamma_analysis_sim)

# Add the executable
add_executable(tas_sim src/tas_sim.cc)
add_executable(ss_sim src/ss_sim.cc)
add_executable(mult_sim src/mult_sim.cc)
add_executable(ss_sim_mult1 src/ss_sim_mult1.cc)
#add_executable(remove_daughter src/remove_daughter.cc)


# Required packages
#find_package(MPI REQUIRED)
#include_directories(${MPI_INCLUDE_PATH})

 
find_package(ROOT REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${ROOT_INCLUDE_DIRS})
target_link_libraries(tas_sim ${ROOT_LIBRARIES}) 
target_link_libraries(ss_sim ${ROOT_LIBRARIES}) 
target_link_libraries(mult_sim ${ROOT_LIBRARIES}) 
target_link_libraries(ss_sim_mult1 ${ROOT_LIBRARIES}) 
#target_link_libraries(remove_daughter ${ROOT_LIBRARIES}) 


# enable warnings
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} `root-config --libs --glibs`")

# enable C++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Enable floating point exception handling 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fnon-call-exceptions -fpermissive")

# generate debugging symbols for release and debug
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")

# the compiler flags
get_directory_property(OUT_VAR LINK_DIRECTORIES) 
message(STATUS "")
message(STATUS "")
message(STATUS "Configuration options:")
message(STATUS "")
message(STATUS "              CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})
message(STATUS "         CMAKE_CXX_LINK_FLAGS: " ${CMAKE_CXX_LINK_FLAGS})
message(STATUS "        CMAKE_CXX_FLAGS_DEBUG: " ${CMAKE_CXX_FLAGS_DEBUG})
message(STATUS "      CMAKE_CXX_FLAGS_RELEASE: " ${CMAKE_CXX_FLAGS_RELEASE})
message(STATUS "              CMAKE_CXX_FLAGS: " ${ROOT_INCLUDE_DIRS})
message(STATUS "             LINK_DIRECTORIES: " ${OUT_VAR})

message(STATUS "")
message(STATUS "")

# Add subdirectories
#add_subdirectory(include)


