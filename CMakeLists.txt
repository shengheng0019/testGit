# cmake_minimum_required(VERSION 3.16)
cmake_minimum_required(VERSION 3.16)
#project(improved)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -DNDEBUG -DIL_STD")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++17")
set(CMAKE_PREFIX_PATH "C:\\CPLEX_Studio1210\\cplex")

set(CMAKE_MOUDLE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
message("${CMAKE_PREFIX_PATH}")
message("1")
set(MY_PROJECT_NAME improved)
project(${MY_PROJECT_NAME})

######## Cplex Dependency

set(CPLEX_ROOT_DIR "C:\\CPLEX_Studio1210")

find_package(cplex)

if(CPLEX_FOUND)
    include_directories(${CPLEX_INCLUDE_DIR})
    link_directories(${CPLEX_LIBRARY})
    message("1")
endif()

INCLUDE_DIRECTORIES("C:\\CPLEX_Studio1210\\cplex\\include")
LINK_DIRECTORIES("C:\\CPLEX_Studio1210\\cplex\\lib")

if(CPLEX_CONCERT_FOUND)
    include_directories(${CPLEX_CONCERT_INCLUDE_DIR})
    link_directories(${CPLEX_CONCERT_LIBRARY})
endif()

INCLUDE_DIRECTORIES("C:\\CPLEX_Studio1210\\concert\\include")
LINK_DIRECTORIES("C:\\CPLEX_Studio1210\\concert\\lib\\x64_windows_msvc14\\stat_mda")

if(CPLEX_ILOCPLEX_FOUND)
    include_directories(${CPLEX_ILOCPLEX_INCLUDE_DIR})
    link_directories(${CPLEX_ILOCPLEX_LIBRARY})
endif()

include_directories("C:\\CPLEX_Studio1210\\cplex\\lib\\x64_windows_msvc14\\stat_mda")
LINK_DIRECTORIES("C:\\CPLEX_Studio1210\\cplex\\lib\\x64_windows_msvc14\\stat_mda")

######## Complier message
set(CMAKE_BUILD_TYPE Release)
message(STATUS "System: ${CMAKE_SYSTEM}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler version: ${CMAKE_CXX_COMPILER_VERSION}")

include_directories("D:/Environment/eigen-3.4.0")
include_directories("D:/Environment/spdlog-1.x/include")
include_directories("D:/Environment/SCIPOptSuite 8.0.3/include")

######## Set output directories.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build//bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build/lib)
link_directories("D:\\Environment\\SCIPOptSuite 8.0.3\\lib")

######## Set C++ 11 standard
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++17" COMPILER_SUPPORTS_CXX17)
if (COMPILER_SUPPORTS_CXX17)
    if (CMAKE_COMPILER_IS_GNUCXX)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++17")
    else ()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
    endif ()
else ()
    message(STATUS "Fatal Error: Your compiler doesn't support c++17, please install it!")
endif ()
file(GLOB SOURCES
        ${PROJECT_SOURCE_DIR}/*.cc
        ${PROJECT_SOURCE_DIR}/*.cpp
        ${PROJECT_SOURCE_DIR}/*.h
        ${PROJECT_SOURCE_DIR}/*.hpp
        "${PROJECT_SOURCE_DIR}/hypergraphpartition/*.hpp"
        "${PROJECT_SOURCE_DIR}/hypergraphpartition/*.cpp"
        )

include_directories(${PROJECT_SOURCE_DIR}/hypergraphpartition)
include_directories(${PROJECT_SOURCE_DIR})

add_executable(${MY_PROJECT_NAME} ${SOURCES})
target_link_libraries(${MY_PROJECT_NAME} cplex12100 ilocplex concert libscip)