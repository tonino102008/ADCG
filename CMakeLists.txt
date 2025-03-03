cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_VERBOSE_MAKEFILE ON)

project(ADCG CXX)

option(BUILD_TESTS "Build test programs." ON)
option(BUILD_IPOPT_TESTS "Build Ipopt test programs." ON)
option(BUILD_LOGGING "Log ADCG and, eventually, the test programs." OFF)
option(BUILD_PROFILER "Profile ADCG and, eventually, the test programs." OFF)
option(BUILD_TIMING "Time ADCG computations." ON)

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    find_package(MKL CONFIG REQUIRED)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")
if (BUILD_IPOPT_TESTS)
    find_package(IPOPT REQUIRED)
endif()

if(BUILD_PROFILER)
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
install(TARGETS ADCG DESTINATION ${INSTALL_DIR})

if(BUILD_TESTS)
    enable_testing()
    set(TEST_DIR "test")
    add_subdirectory(${TEST_DIR})
endif()