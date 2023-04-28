cmake_minimum_required(VERSION 3.9)
project(RangeKNN)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# MESSAGE(STATUS ${OpenMP_CXX})

set(CMAKE_CXX_STANDARD 11)

message(${CMAKE_SYSTEM})


find_package(Boost REQUIRED COMPONENTS timer chrono system program_options)


message(STATUS ${Boost_LIBRARIES})


if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
    SET (LD_LIBRARY_PATH ${BOOST_LIBRARYDIR})
endif()

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src)

set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -O3 ")

option(NOSIMD "build hnsw on SIMD" ON)
IF(NOSIMD)
    add_definitions(-DNO_MANUAL_VECTORIZATION)
    add_definitions(-DKGRAPH_NO_VECTORIZE)
ENDIF(NOSIMD)

add_subdirectory(src/utils)
add_subdirectory(include/kgraph)
add_subdirectory(tests)

