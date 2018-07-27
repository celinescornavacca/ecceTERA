# CMake script for Bio++ and Treerecs Library
# Authors:
#   Nicolas COMTE
# Created: 27/07/2017

cmake_minimum_required(VERSION 3.0.1)

# Add Bpp sources
add_subdirectory(bpp-core/src)
add_subdirectory(bpp-seq/src)
add_subdirectory(bpp-phyl/src)

# Export variables BPP_INCLUDE_DIR and BPP_LIBRARY in parent scope
set(BPP_INCLUDE_DIR ${BPP_INCLUDE_DIR} PARENT_SCOPE)
set(BPP_LIBRARY ${BPP_LIBRARY} PARENT_SCOPE)