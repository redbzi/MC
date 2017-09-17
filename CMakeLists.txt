cmake_minimum_required(VERSION 3.6)
project(SM)

SET(CMAKE_CXX_STANDARD 11)

SET(SOURCE_FILES
        main.cpp
        basic.cpp
        analytic.cpp
        heston.cpp
        plot.cpp
        print.cpp)

SET(LIBRARIES
        libboost_system.dylib
        libboost_iostreams.dylib
        libboost_filesystem.dylib)

add_definitions("-Wall -O3")

ADD_EXECUTABLE(SM_Heston ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(SM_Heston ${LIBRARIES})