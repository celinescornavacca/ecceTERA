
include_directories(${BPP_INCLUDE_DIR})

set (BPP_PHYL_CPP_FILES
Bpp/Phyl/BipartitionList.cpp
Bpp/Phyl/NNITopologySearch.cpp
Bpp/Phyl/BipartitionTools.cpp
Bpp/Phyl/OptimizationTools.cpp
Bpp/Phyl/Likelihood/PseudoNewtonOptimizer.cpp
Bpp/Phyl/Likelihood/DRHomogeneousTreeLikelihood.cpp
Bpp/Phyl/Likelihood/DRASDRTreeLikelihoodData.cpp
Bpp/Phyl/Likelihood/AbstractHomogeneousTreeLikelihood.cpp
Bpp/Phyl/Likelihood/AbstractDiscreteRatesAcrossSitesTreeLikelihood.cpp
Bpp/Phyl/Likelihood/AbstractTreeLikelihood.cpp
Bpp/Phyl/Likelihood/NNIHomogeneousTreeLikelihood.cpp
Bpp/Phyl/Likelihood/GlobalClockTreeLikelihoodFunctionWrapper.cpp
Bpp/Phyl/Io/Newick.cpp
Bpp/Phyl/Node.cpp
Bpp/Phyl/TreeExceptions.cpp
Bpp/Phyl/TreeTemplateTools.cpp
Bpp/Phyl/TreeTools.cpp
Bpp/Phyl/Distance/BioNJ.cpp
Bpp/Phyl/Distance/DistanceEstimation.cpp
Bpp/Phyl/Distance/AbstractAgglomerativeDistanceMethod.cpp
Bpp/Phyl/Distance/NeighborJoining.cpp
Bpp/Phyl/PatternTools.cpp
Bpp/Phyl/Parsimony/DRTreeParsimonyScore.cpp
Bpp/Phyl/Parsimony/DRTreeParsimonyData.cpp
Bpp/Phyl/Parsimony/AbstractTreeParsimonyScore.cpp
Bpp/Phyl/Model/StateMap.cpp
Bpp/Phyl/Model/AbstractSubstitutionModel.cpp
Bpp/Phyl/Model/Nucleotide/JCnuc.cpp
Bpp/Phyl/SitePatterns.cpp)

# Build the static lib
add_library (bpp-phyl STATIC ${BPP_PHYL_CPP_FILES})

# Create library
target_link_libraries (bpp-phyl ${BPP_LIBRARY})

# Update variables and set for parent scope
set(BPP_INCLUDE_DIR ${BPP_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/ext/bpp/bpp-phyl/src PARENT_SCOPE)
set(BPP_LIBRARY ${BPP_LIBRARY} bpp-phyl PARENT_SCOPE)