cmake_minimum_required(VERSION 3.15.0)
project(SEA VERSION 0.1.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 14)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Do not use shared lib.")


include(FetchContent)
include(cmake/Eigen3.cmake)
include(cmake/cli11.cmake)
include(cmake/nlopt.cmake)


aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC_LIST)

add_library(SEA STATIC ${SRC_LIST})
target_include_directories(SEA PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(SEA PUBLIC Eigen3::Eigen)

# OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(SEA PUBLIC OpenMP::OpenMP_CXX)
endif()


add_executable(SEA_QN apps/SEA_QN.cpp)
target_link_libraries(SEA_QN PUBLIC SEA CLI11::CLI11 nlopt::nlopt)
