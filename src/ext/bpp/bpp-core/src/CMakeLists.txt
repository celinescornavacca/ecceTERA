
# File list
set (BPP_CORE_CPP_FILES
Bpp/App/ApplicationTools.cpp
Bpp/Numeric/Function/AbstractOptimizer.cpp
Bpp/Numeric/Function/BrentOneDimension.cpp
Bpp/Numeric/Function/ConstantDistribution.cpp
Bpp/Numeric/Function/DownhillSimplexMethod.cpp
Bpp/Numeric/Function/AbstractDiscreteDistribution.cpp
Bpp/Numeric/Function/OneDimensionOptimizationTools.cpp
Bpp/Numeric/Function/DirectionFunction.cpp
Bpp/Numeric/Function/TwoPointsNumericalDerivative.cpp
Bpp/Numeric/Function/ThreePointsNumericalDerivative.cpp
Bpp/Numeric/Function/NewtonBacktrackOneDimension.cpp
Bpp/Numeric/Function/ConjugateGradientMultiDimensions.cpp
Bpp/Numeric/Function/ReparametrizationFunctionWrapper.cpp
Bpp/BppString.cpp
Bpp/Numeric/Function/BfgsMultiDimensions.cpp
Bpp/Exceptions.cpp
Bpp/Io/FileTools.cpp
Bpp/Numeric/AbstractParameterAliasable.cpp
Bpp/Numeric/Function/OptimizationStopCondition.cpp
Bpp/Numeric/Function/SimpleMultiDimensions.cpp
Bpp/Numeric/Function/MetaOptimizer.cpp
Bpp/Numeric/AbstractParametrizable.cpp
Bpp/Numeric/AdaptiveKernelDensityEstimation.cpp
Bpp/Numeric/NumTools.cpp
Bpp/Numeric/Parameter.cpp
Bpp/Numeric/AutoParameter.cpp
Bpp/Numeric/ParameterExceptions.cpp
Bpp/Numeric/ParameterList.cpp
Bpp/Numeric/Random/RandomTools.cpp
Bpp/Numeric/Random/Uniform01K.cpp
Bpp/Numeric/VectorTools.cpp
Bpp/Text/NestedStringTokenizer.cpp
Bpp/Text/StringTokenizer.cpp
Bpp/Text/TextTools.cpp
Bpp/Utils/AttributesTools.cpp
        )

# Build the static lib
add_library (bpp-core STATIC ${BPP_CORE_CPP_FILES})

# Add current dir into BPP_INCLUDE_DIR variable and library into BPP_LIBRARY
# then set for the parent_scope
set(BPP_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/ext/bpp/bpp-core/src PARENT_SCOPE)
set(BPP_LIBRARY bpp-core PARENT_SCOPE)