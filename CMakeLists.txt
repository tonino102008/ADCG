cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_VERBOSE_MAKEFILE ON)

project(indOptLap LANGUAGES C CXX)

option(BUILD_TESTS "Build test programs." ON)
option(BUILD_LOGGING "Log indOptLap and, eventually, the test programs." OFF)
option(BUILD_PROFILER "Profile indOptLap and, eventually, the test programs." OFF)
option(BUILD_TIMING "Time indOptLap linear system solver." OFF)

find_package(MKL CONFIG REQUIRED)
set(CMAKE_C_COMPILER icx)
set(CMAKE_CXX_COMPILER icpx)


if(BUILD_PROFILER)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif()

if(BUILD_LOGGING)
    add_definitions(-DLOGGING)
endif()

if(BUILD_TIMING)
    add_definitions(-DTIMING)
endif()

set(SRC_DIR "src")
add_subdirectory(${SRC_DIR})

set(INSTALL_DIR "bin")
install(TARGETS indOptLap DESTINATION ${INSTALL_DIR})

if(BUILD_TESTS)
    enable_testing()
    set(TEST_DIR "test")
    add_subdirectory(${TEST_DIR})
endif()