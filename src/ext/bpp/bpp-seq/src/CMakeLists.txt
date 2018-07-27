
include_directories(${BPP_INCLUDE_DIR})

# File list
SET(BPP_SEQ_CPP_FILES
Bpp/Seq/Alphabet/AlphabetExceptions.cpp
Bpp/Seq/Alphabet/AlphabetTools.cpp
Bpp/Seq/Alphabet/DNA.cpp
Bpp/Seq/Alphabet/RNA.cpp
Bpp/Seq/Alphabet/WordAlphabet.cpp
Bpp/Seq/Alphabet/LetterAlphabet.cpp
Bpp/Seq/Alphabet/ProteicAlphabet.cpp
Bpp/Seq/Alphabet/DefaultAlphabet.cpp
Bpp/Seq/Alphabet/AbstractAlphabet.cpp
Bpp/Seq/Alphabet/CodonAlphabet.cpp
Bpp/Seq/Container/AbstractSequenceContainer.cpp
Bpp/Seq/Container/AlignedSequenceContainer.cpp
Bpp/Seq/Container/VectorSequenceContainer.cpp
Bpp/Seq/Container/VectorSiteContainer.cpp
Bpp/Seq/Container/SequenceContainerTools.cpp
Bpp/Seq/DistanceMatrix.cpp
Bpp/Seq/SiteTools.cpp
Bpp/Seq/Sequence.cpp
Bpp/Seq/SequenceExceptions.cpp
Bpp/Seq/SymbolListTools.cpp
Bpp/Seq/Site.cpp
Bpp/Seq/SiteExceptions.cpp
Bpp/Seq/StringSequenceTools.cpp
Bpp/Seq/SymbolList.cpp
        )

# Build the static lib
add_library (bpp-seq STATIC ${BPP_SEQ_CPP_FILES})

# Link bpp-core
target_link_libraries (bpp-seq ${BPP_LIBRARY})

# Add bpp-seq sources dir into BPP_INCLUDE_DIR and lib into BPP_LIBRARY, set for PARENT SCOPE
set(BPP_INCLUDE_DIR ${BPP_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/ext/bpp/bpp-seq/src PARENT_SCOPE)
set(BPP_LIBRARY ${BPP_LIBRARY} bpp-seq PARENT_SCOPE)