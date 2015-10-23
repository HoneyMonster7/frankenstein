cmake_minimum_required (VERSION 2.6)

project (Reaction)

find_package(Boost REQUIRED)

add_executable(reaction main.cpp reaction.cpp)

include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(reaction ${Boost_LIBRARIES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")
add_definitions("-DDATA_PATH=\"${CMAKE_SOURCE_DIR}/../data/\"")
