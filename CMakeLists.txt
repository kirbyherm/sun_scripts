
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)

#set(CMAKE_VERBOSE_MAKEFILE ON) 

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules/")
message(STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH}) 
project(gamma_analysis)

# Add the executable
add_executable(tas src/tas.cc)
add_executable(tas_ss_2d src/tas_ss_2d.cc)
add_executable(ss src/ss.cc)
add_executable(mult src/mult.cc)
add_executable(remove_daughter src/remove_daughter.cc)
add_executable(write_to_csv src/write_to_csv.cc)


# Required packages
#find_package(MPI REQUIRED)
#include_directories(${MPI_INCLUDE_PATH})

 
find_package(ROOT REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${ROOT_INCLUDE_DIRS})
target_link_libraries(tas ${ROOT_LIBRARIES}) 
target_link_libraries(tas_ss_2d ${ROOT_LIBRARIES}) 
target_link_libraries(ss ${ROOT_LIBRARIES}) 
target_link_libraries(mult ${ROOT_LIBRARIES}) 
target_link_libraries(remove_daughter ${ROOT_LIBRARIES}) 
target_link_libraries(write_to_csv ${ROOT_LIBRARIES}) 


# enable warnings
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")

# enable C++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Enable floating point exception handling 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fnon-call-exceptions")

# generate debugging symbols for release and debug
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

# enable root libs
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

install(TARGETS tas DESTINATION bin)
install(TARGETS tas_ss_2d DESTINATION bin)
install(TARGETS ss DESTINATION bin)
install(TARGETS mult DESTINATION bin)
install(TARGETS remove_daughter DESTINATION bin)
install(TARGETS write_to_csv DESTINATION bin)
