# Minimum version of CMake required
cmake_minimum_required(VERSION 3.10)

# Project name and version
project(pimavilo_mapper VERSION 1.0.0)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#Add subdirectory bioparser
add_subdirectory(bioparser)

# Add the executable
add_executable(msa src/msa.cpp)

# Link bioparser and pimavilo to the executable
target_link_libraries(msa
    PRIVATE bioparser
)