# Main code
add_library(mma_mma MMASolver.cpp)
target_include_directories(mma_mma PUBLIC ..)

# Let's get a little bit paranoid
include(SetWarnings)
target_compile_options(mma_mma PRIVATE ${ALL_WARNINGS})

# Use C++11
include(CXXFeatures)
target_compile_features(mma_mma PUBLIC ${CXX11_FEATURES})

# Alias target
add_library(mma::mma ALIAS mma_mma)
