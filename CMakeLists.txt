CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(ecceTERA CXX)

# Change project root and source vars
set(PROJECT_ROOT_DIR ${PROJECT_SOURCE_DIR})

set(PROJECT_SOURCE_DIR ${PROJECT_ROOT_DIR}/src)

# Tell CMake where to look for custom modules
set(CMAKE_MODULE_PATH ${PROJECT_ROOT_DIR}/cmake)

# Targets
add_subdirectory(src) #add sources (Treerecs + misc + Bpp + Boost)

