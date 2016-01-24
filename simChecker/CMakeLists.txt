cmake_minimum_required (VERSION 2.6)

project (simMatrix)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../cmake/")

#find_package(Boost REQUIRED)
#find_package(GLPK REQUIRED)

add_executable(simMatrix main.cpp )
#add_executable(reaction main.cpp reaction.cpp cell cell.cpp)

#include_directories(${Boost_INCLUDE_DIRS})
#target_link_libraries(reaction ${Boost_LIBRARIES})


#include_directories(${GLPK_INCLUDE_DIRS})
#target_link_libraries(reaction ${GLPK_LIBRARIES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")
#add_definitions("-DDATA_PATH=\"${CMAKE_SOURCE_DIR}/../data/\"")
