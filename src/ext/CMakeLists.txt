
# Add _USE_MATH_DEFINES for cmath.
add_definitions(-D_USE_MATH_DEFINES)

# Add Bio++ and Boost sources and create a libraries
add_subdirectory(boost)
add_subdirectory(bpp)


# Set variables into cmake parent scope.
# Bio++ core/seq/phyl
set(BPP_INCLUDE_DIR ${BPP_INCLUDE_DIR} PARENT_SCOPE)
set(BPP_LIBRARY ${BPP_LIBRARY} PARENT_SCOPE)

# Boost
set(BOOST_INCLUDE_DIR ${BOOST_INCLUDE_DIR} PARENT_SCOPE)

