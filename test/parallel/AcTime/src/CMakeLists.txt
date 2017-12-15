cmake_minimum_required (VERSION 2.6)
project (ParAcTest)

find_package(Threads REQUIRED)

# include the file where I define a couple of compiler flags, including the long
# lists of warning flags for the compiler
include(./CMakeFlagDefs.txt)

# add our default g++ flags
#add_definitions(-m64 -std=c++14 -D_FILE_OFFSET_BITS=64 -lpthread -I/usr/local/include/root)
add_definitions(-m64 -march=native -mtune=native -std=c++11 -D_FILE_OFFSET_BITS=64)

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

# set things up so that we pull in the MCMCpp headers
include_directories("../../../../MCMCpp")

# Set the executable and what it depends on
add_executable(acTest ${SRCS} ${HDRS})
target_link_libraries(acTest ${CMAKE_THREAD_LIBS_INIT})
