SET(CMAKE_INCLUDE_CURRENT_DIR ON)
ADD_EXECUTABLE( ecceTERA DTLMatrixNetwork.cpp ecceTERA.cpp DTLMatrix.cpp CladesAndTripartitions.cpp MyMatrix.cpp MySpeciesTree.cpp MyGeneTree.cpp DTLGraph.cpp MyMatrixV.cpp MyMatrixT.cpp DTLMatrixTriplets.cpp DTLMatrixRecalc.cpp MyNetwork.cpp NetAlg.cpp  )
target_link_libraries( ecceTERA ${LIBS} )
install(TARGETS ecceTERA DESTINATION bin )
