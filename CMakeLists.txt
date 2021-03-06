cmake_minimum_required(VERSION 3.16)
project(f4d)

set(CMAKE_CXX_STANDARD 17)

set(WARNING_FLAGS "-Wall -Wextra -Wpedantic -Wuninitialized -Wshadow")
set(SANITIZERS_FLAGS "-fno-omit-frame-pointer -fsanitize=address -fsanitize-address-use-after-scope -fsanitize=undefined")
set(OPT_FLAGS "-O3 -march=native -ffat-lto-objects -flto")


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WARNING_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${SANITIZERS_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OPT_FLAGS}")

set(CMAKE_VERBOSE_MAKEFILE ON)

set(SOURCE main.cpp
        greedy_bpp.hpp 
        routemin.hpp  
        RuinAndRecreate.hpp 
        LinKernighan.hpp  
        greedy_bpp.cpp)

# Find CPLEX
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(CPLEX)
if (CPLEX_FOUND)
    include_directories(${CPLEX_INCLUDE_DIRS})
    message(${CPLEX_INCLUDE_DIRS})
    set(LIBRARIES ${LIBRARIES} ${CPLEX_LIBRARIES} pthread dl m)
    add_definitions(-DIL_STD)
endif(CPLEX_FOUND)

add_executable(f4d ${SOURCE})


target_link_libraries(f4d PUBLIC ${LIBRARIES} stdc++fs)