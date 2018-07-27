cmake_minimum_required(VERSION 3.0.1)
add_subdirectory(src)

# Export variables BPP_INCLUDE_DIR and BPP_LIBRARY in parent scope
set(BPP_INCLUDE_DIR ${BPP_INCLUDE_DIR} PARENT_SCOPE)
set(BPP_LIBRARY ${BPP_LIBRARY} PARENT_SCOPE)