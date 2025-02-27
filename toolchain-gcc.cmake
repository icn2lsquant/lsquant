set(CMAKE_CXX_STANDARD 11) # tODO move up to a general cmake config for all sub projects ?

set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_COMPILER_VENDOR "intel")

set(CMAKE_C_COMPILER gcc)
set(CMAKE_C_FLAGS " -O3 --std=c++11 -lgsl -Wall -Wextra -fopenmp " CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELEASE " ")
set(CMAKE_C_FLAGS_DEBUG "-g -fopenmp")

set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "-O3 --std=c++11 -lgsl -Wall -Wextra -fopenmp " CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "")
set(CMAKE_CXX_FLAGS_DEBUG "-g -fopenmp")

set(INTEL_MKL "FALSE")

set(EIGEN "TRUE")
set(BLAS_LIB "-lpthread -lfftw3 -lgsl -lm -ldl")

