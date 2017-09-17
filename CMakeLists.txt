cmake_minimum_required(VERSION 3.8)
project(MC-Heston)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_BUILD_TYPE release)

set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})

include(FindBoost)

find_package(
        Boost

        1.63.0

        REQUIRED
)

link_directories(${Boost_LIBRARY_DIRS})

include_directories(${Boost_INCLUDE_DIRS})

set(SOURCE_FILES
        main
        basic
        analytic
        heston
        plot
        print
        )

set(LIBRARIES
        boost_system
        boost_iostreams
        boost_filesystem
        )

add_definitions("-Wall -O3")

add_executable(run ${SOURCE_FILES})
target_link_libraries(run ${LIBRARIES})
