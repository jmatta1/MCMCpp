#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
cmake_minimum_required (VERSION 3.8)
project (SeqSkewedGaussMetHast)

set(MCMCPP_DIR ../../../../../MCMCpp)
set(COMMON_DIR ../../Common)

# include the file where I define a couple of compiler flags, including the long
# lists of warning flags for the compiler
include(${MCMCPP_DIR}/MCMCppFlags.cmake)
include(./CMakeFlagDefs.txt)

# add our default g++ flags
#add_definitions(-m64 -std=c++14 -D_FILE_OFFSET_BITS=64 -lpthread -I/usr/local/include/root)
#add_definitions(-std=c++14 -stdlib=libc++ -m64 -march=native -mtune=native -D_FILE_OFFSET_BITS=64)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -march=native -mtune=native")

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

# Set the source file
set(SRC "main.cpp")

# Set the objective function header file
set(HDR "${COMMON_DIR}/SkewedGaussian.h")

# set things up so that we pull in the MCMCpp headers into the Qt editor
file(GLOB_RECURSE MCMCPP_HDRS ${MCMCPP_DIR}/*.h)

# set things up so that files include as if they were at the root of the dir
# tree unless they are including a file within their own folder
include_directories(${COMMON_DIR})
include_directories(${MCMCPP_DIR})

# Set the executable and what it depends on
add_executable(metHast ${SRC} ${HDR} ${MCMCPP_HDRS})

