cmake_minimum_required(VERSION 3.8)
project(MC)

set(CMAKE_CXX_STANDARD 11)

include(FindBoost)

find_package(
        Boost

        1.63.0

        REQUIRED
)

link_directories(${Boost_LIBRARY_DIRS})

include_directories(${Boost_INCLUDE_DIRS})

set(SOURCE
        gbm
        heston
        )

set(LIBRARIES
        boost_system
        boost_iostreams
        boost_filesystem
        )

add_definitions("-Wall -O2")

add_library(mc SHARED ${SOURCE})
target_link_libraries(mc ${LIBRARIES})
