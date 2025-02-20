find_path(IPOPT_INCLUDE_DIR NAMES IpNLP.hpp
    PATHS   "/usr/include/coin"
            "/usr/include/coin-or"
            "/usr/local/include/coin"
            "/usr/local/include/coin-or"
)

find_path(IPOPT_MUMPS_INCLUDE_DIR NAMES dmumps_c.h
    PATHS   "/usr/include/coin"
            "/usr/include/coin-or"
            "/usr/local/include/coin"
            "/usr/local/include/coin-or"
            "/usr/local/include/coin/mumps"
            "/usr/local/include/coin-or/mumps"
)

find_library(IPOPT_LIBRARY ipopt
    PATHS
    "/usr/lib"
    "/usr/local/lib"
)

find_path(IPOPT_MUMPS_LIBRARY_DIR 
    NAMES libcoinmumps.so lapack -lblas -lm  -ldl
    PATHS 
    "/usr/lib"
    "/usr/local/lib"
)

find_library(IPOPT_MUMPS_LIBRARY coinmumps lapack blas
    PATHS 
    "/usr/lib"
    "/usr/local/lib"
)

set(IPOPT_LIBRARY_DIRS "${IPOPT_LIBRARY_DIR}")
set(IPOPT_INCLUDE_DIRS "${IPOPT_INCLUDE_DIR}")
set(IPOPT_LIBRARIES "${IPOPT_LIBRARY}")

if(IPOPT_MUMPS_INCLUDE_DIR)
    list(APPEND IPOPT_INCLUDE_DIR "${IPOPT_MUMPS_INCLUDE_DIR}")   
endif(IPOPT_MUMPS_INCLUDE_DIR)

if(IPOPT_MUMPS_LIBRARY_DIR)
    message ("IPOPT_MUMPS_LIBRARY_DIR found at ${IPOPT_MUMPS_LIBRARY_DIR}")
    list(APPEND IPOPT_MUMPS_LIBRARY "${IPOPT_MUMPS_LIBRARY_DIR}")
endif(IPOPT_MUMPS_LIBRARY_DIR)

if(IPOPT_MUMPS_LIBRARY)
    list(APPEND IPOPT_LIBRARIES "${IPOPT_MUMPS_LIBRARY}")   
endif(IPOPT_MUMPS_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(IPOPT DEFAULT_MSG IPOPT_LIBRARY IPOPT_INCLUDE_DIR)

mark_as_advanced(IPOPT_INCLUDE_DIR IPOPT_LIBRARY)
