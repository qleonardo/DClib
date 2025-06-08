cmake_minimum_required (VERSION 2.6)

# Compilers
set (LIBRARY_OUTPUT_PATH ../lib/${CMAKE_BUILD_TYPE})
set (CMAKE_CXX_COMPILER icpc)

# Project
project (TDTree-lib)

add_definitions (-D${SHARED})
if (${FORKJOIN})
    add_definitions (-DFORKJOIN)
endif (${FORKJOIN})

if (${OMP})
    set (flags "${flags} -qopenmp")
    add_definitions (-DOMP)
endif (${OMP})

set (CMAKE_CXX_FLAGS "${flags} -O3 -std=c++11")

# Headers
include_directories (../include)
include_directories (../src/headers)

# Sources
file (GLOB_RECURSE source_list ../src/* ../include/*)

# Library
set (lib DC)

add_library (
	${lib}
    SHARED
	${source_list}
)
