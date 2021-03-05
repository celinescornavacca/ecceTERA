SET(CMAKE_INCLUDE_CURRENT_DIR ON)



## Add targets
# Add Bio++ and Boost sources and create a libraries
add_subdirectory(ext)



ADD_EXECUTABLE( ecceTERA ecceTERA.cpp DTLMatrix.cpp CladesAndTripartitions.cpp MyMatrix.cpp MySpeciesTree.cpp MyGeneTree.cpp DTLGraph.cpp )

target_include_directories(${PROJECT_NAME} PRIVATE ${BPP_INCLUDE_DIR} ${BOOST_INCLUDE_DIR} ${SIMPLESVG_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} ${BPP_LIBRARY})
set_target_properties(${PROJECT_NAME}
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        )


# Add install target
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY true)
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)


## Export variables into parent scope
# Bio++ core/seq/phyl
set(BPP_INCLUDE_DIR ${BPP_INCLUDE_DIR} PARENT_SCOPE)
set(BPP_LIBRARY ${BPP_LIBRARY} PARENT_SCOPE)

# Boost
set(BOOST_INCLUDE_DIR ${BOOST_INCLUDE_DIR} PARENT_SCOPE)
