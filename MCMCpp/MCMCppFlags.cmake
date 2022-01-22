# Handle getting the C++ standard to C++11
if (CMAKE_VERSION VERSION_LESS "3.1")
    # if the cmake version  is less than 3.1 we have to do this the hard way
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lc++")
    endif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
else ()
    # Otherwise we can do this the easy way
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif ()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_FILE_OFFSET_BITS=64")

if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DOSX_IS_ANNOYING")
endif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")


if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGCC_SUPPORTS_FANCYNESS")
endif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
