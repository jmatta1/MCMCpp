cmake_minimum_required (VERSION 2.6)
project (SeqSkewedGaussMetHast)

# include the file where I define a couple of compiler flags, including the long
# lists of warning flags for the compiler
include(./CMakeFlagDefs.txt)
include("../../../../../MCMCpp/CMakeMCMCppFlags.txt")

# add our default g++ flags
#add_definitions(-m64 -std=c++14 -D_FILE_OFFSET_BITS=64 -lpthread -I/usr/local/include/root)
#add_definitions(-std=c++14 -stdlib=libc++ -m64 -march=native -mtune=native -D_FILE_OFFSET_BITS=64)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -march=native -mtune=native -D_FILE_OFFSET_BITS=64")

# a little hack to make Qt Creator a little happier when it runs cmake to get
# the project files
# check for each supported build type and add the appropriate g++ flags
if(BUILD_TYPE STREQUAL "Release")
    add_definitions(${GCC_OPT_FLAG})
elseif(BUILD_TYPE STREQUAL "Debug")
    add_definitions(${GCC_DEB_FLAG})
elseif(BUILD_TYPE STREQUAL "OptDebug")
    add_definitions(${GCC_DEB_FLAG})
    add_definitions(${GCC_OPT_FLAG})
elseif(BUILD_TYPE STREQUAL "Warn")
    add_definitions(${WARN_FLAGS})
elseif(BUILD_TYPE STREQUAL "Warn_Opt")
    add_definitions(${WARN_FLAGS})
    add_definitions(${WARN_OPT_FLAGS})
endif(BUILD_TYPE STREQUAL "Release")

# Grab all the source files
file(GLOB_RECURSE SRCS *.cpp)

# Grab all the header files
file(GLOB_RECURSE HDRS *.h)

# set things up so that files include as if they were at the root of the dir
# tree unless they are including a file within their own folder
include_directories(".")
include_directories("../../../../../MCMCpp")

# set things up so that we pull in the MCMCpp headers into the Qt editor
#include_directories("../../../../../MCMCpp")
file(GLOB_RECURSE HDRS ../../../../../MCMCpp/*.h)

# Set the executable and what it depends on
add_executable(metHast ${SRCS} ${HDRS})

