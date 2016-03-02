cmake_minimum_required (VERSION 2.6)


project (Reaction)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../cmake/")

find_package(Boost REQUIRED)
find_package(GLPK REQUIRED)

SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
SET(BUILD_SHARED_LIBRARIES OFF)
SET(CMAKE_EXE_LINKER_FLAGS "-static")
#SET(CMAKE_SHARED_LINKER_FLAGS "-lltdl -ldl -lrt")

add_executable(reaction main.cpp reaction.cpp substrate.cpp cell cell.cpp)
#add_executable(reaction main.cpp reaction.cpp cell cell.cpp)

include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(reaction ${Boost_LIBRARIES})


include_directories(${GLPK_INCLUDE_DIRS})
target_link_libraries(reaction -lltdl  ${GLPK_LIBRARIES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")
add_definitions("-DDATA_PATH=\"${CMAKE_SOURCE_DIR}/../data/\"")

set(CMAKE_BUILD_TYPE Debug)

