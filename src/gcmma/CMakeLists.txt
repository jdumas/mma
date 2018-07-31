# Main code
add_library(mma_gcmma GCMMASolver.cpp)
target_include_directories(mma_gcmma PUBLIC ..)

# Let's get a little bit paranoid
include(SetWarnings)
target_compile_options(mma_gcmma PRIVATE ${ALL_WARNINGS})

# Use C++11
include(CXXFeatures)
target_compile_features(mma_gcmma PUBLIC ${CXX11_FEATURES})

# Alias target
add_library(mma::gcmma ALIAS mma_gcmma)
