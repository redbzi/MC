cmake_minimum_required(VERSION 3.8)
project(MC)

set(CMAKE_CXX_STANDARD 11)

link_directories(../cmake-build-debug/src/)

set(SOURCE
        interface
        plot
        print
        )

set(LIBRARIES
        mc
        )

add_definitions("-Wall -O2")

add_executable(run ${SOURCE})
target_link_libraries(run ${LIBRARIES})
