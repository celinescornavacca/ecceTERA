CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(tera CXX)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)


SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# --- the IMPROVED_FIND_LIBRARY macro (taken from the Bio++ CMakeLists.txt) ---
macro(IMPROVED_FIND_LIBRARY OUTPUT_LIBS lib_name include_to_find)
    FIND_PATH(${lib_name}_INCLUDE_DIR ${include_to_find} PATH_SUFFIXES bpp-phyl-2.2.0/src bpp-seq-2.2.0/src bpp-core-2.2.0/src )
    SET(${lib_name}_NAMES ${lib_name} ${lib_name}lib ${lib_name}dll)
    FIND_LIBRARY(${lib_name}_LIBRARY NAMES ${${lib_name}_NAMES} PATH_SUFFIXES bpp-phyl-2.2.0/src bpp-seq-2.2.0/src bpp-core-2.2.0/src )

    IF(${lib_name}_LIBRARY)
        MESSAGE("-- Library ${lib_name} found here:")
        MESSAGE("    ${${lib_name}_INCLUDE_DIR}/")
        MESSAGE("    ${${lib_name}_LIBRARY}")
    ELSE(${lib_name}_LIBRARY)
        MESSAGE(FATAL_ERROR "${lib_name} required but not found.")
    ENDIF(${lib_name}_LIBRARY)

    #add the dependency:
    INCLUDE_DIRECTORIES(${${lib_name}_INCLUDE_DIR})
    SET(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_LIBRARY})
ENDMACRO(IMPROVED_FIND_LIBRARY)

#--- Libraries & includes ---
include_directories("include/")

# uncomment for static compilation
#set(CMAKE_FIND_LIBRARY_SUFFIXES ".a ${CMAKE_FIND_LIBRARY_SUFFIXES}") 

IMPROVED_FIND_LIBRARY(LIBS bpp-phyl Bpp/Phyl/Tree.h)
IMPROVED_FIND_LIBRARY(LIBS bpp-seq Bpp/Seq/Alphabet/Alphabet.h)
IMPROVED_FIND_LIBRARY(LIBS bpp-core Bpp/Clonable.h)

#Boost libraries:
#find_package(Boost COMPONENTS graph REQUIRED)
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

include_directories("include")

# Source Directory
add_subdirectory(src)

