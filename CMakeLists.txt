cmake_minimum_required(VERSION 2.8)
project(arap_morph)

file(GLOB headers "${CMAKE_SOURCE_DIR}/src/*.h")
file(GLOB sources "${CMAKE_SOURCE_DIR}/src/*.cpp")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(vdb openvdb tbb Half)
set(openmesh OpenMeshCore)
set(igl igl)
set(lemon emon)

set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE " ${CMAKE_CXX_FLAGS_RELEASE} -Wreturn-type -fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -Wno-sign-compare")
#message("CMAKE_CXX_FLAGS_RELEASE is " ${CMAKE_CXX_FLAGS_RELEASE})

include_directories(/usr/include/eigen3/ /usr/local/OpenVDB/include/)
link_directories(/usr/local/OpenVDB/lib/ /usr/local/lib/igl/)

add_executable(morph ${headers} ${sources})
target_link_libraries(morph ${vdb} ${openmesh} ${igl} ${lemon})

