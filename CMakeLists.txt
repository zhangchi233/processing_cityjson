cmake_minimum_required (VERSION 3.14)
project (calcu)

set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_FLAGS "-Wall -std=c++14")
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-O2" )

# CGAL
find_package( CGAL QUIET COMPONENTS )
if ( CGAL_FOUND )
    include( ${CGAL_USE_FILE} )
    message(STATUS ${CGAL_LIBRARIES})
    message(STATUS ${CGAL_3RD_PARTY_LIBRARIES})
else()
    message(SEND_ERROR "CGAL library is required")
    return()
endif()

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include(CGAL_Eigen3_support)

include_directories( ${CMAKE_SOURCE_DIR}/include/ )

FILE(GLOB SRC_FILES src/*.cpp)
add_executable(calcu ${SRC_FILES})

target_link_libraries(${PROJECT_NAME} CGAL::CGAL CGAL::Eigen3_support)

