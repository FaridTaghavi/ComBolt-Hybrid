# CMakeLists.txt

cmake_minimum_required(VERSION 3.10)
project(ComBolt-ITA LANGUAGES CXX)

# Compiler settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

# ----- ROOT -----
set(USE_ROOT TRUE)
if (USE_ROOT)
    message(STATUS "Looking for ROOT")
    find_package(ROOT 6.21 CONFIG REQUIRED)
    message(STATUS "ROOT ${ROOT_VERSION} found at ${ROOT_DIR}.")
    include(${ROOT_USE_FILE})
endif ()

# ----- Boost -----
message(STATUS "Looking for BOOST library")
message(STATUS "Looking for BOOST library in ${BOOST_ROOT}")

# Use system Boost ≥1.75 (your node has 1.75). If you later install 1.87 in $HOME,
# you can configure with:
#   -DBOOST_ROOT=$HOME/opt/boost-1.87.0 -DBoost_NO_SYSTEM_PATHS=ON -DCMAKE_PREFIX_PATH=$HOME/opt/boost-1.87.0
find_package(Boost 1.75 REQUIRED COMPONENTS filesystem program_options system)

message(STATUS "Boost include directories: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost library directories: ${Boost_LIBRARY_DIRS}")

# ----- HDF5 -----
message(STATUS "Looking for HDF5 library")
find_package(HDF5 COMPONENTS C CXX REQUIRED)
message(STATUS "HDF5 include directory: ${HDF5_INCLUDE_DIRS}")
message(STATUS "HDF5 libraries: ${HDF5_LIBRARIES}")

# ----- GSL -----
message(STATUS "Looking for GSL library")
find_package(GSL REQUIRED)
message(STATUS "GSL include dir found at ${GSL_INCLUDE_DIR}")

# ----- OpenMP -----
message(STATUS "Looking for OpenMP library")
find_package(OpenMP REQUIRED)
message(STATUS "OpenMP found")

# Include paths (use target_include_directories for new targets where possible)
include_directories(
  ${CMAKE_SOURCE_DIR}/src
  SYSTEM
  ${Boost_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${GSL_INCLUDE_DIR}
)

# Your project sources
add_subdirectory(src)

# Main executable
add_executable(ComBolt-ITA ComBolt-ITA.cpp)

# Link libraries
target_link_libraries(ComBolt-ITA
  PRIVATE
    src
    ${ROOT_LIBRARIES}
    Boost::filesystem
    Boost::program_options
    Boost::system
    ${HDF5_LIBRARIES}
    GSL::gsl
    GSL::gslcblas
    OpenMP::OpenMP_CXX
)

# Note:
# - stdc++fs is not needed on GCC 11 with C++17; add it back ONLY if your linker complains:
#     target_link_libraries(ComBolt-ITA PRIVATE stdc++fs)
# - If you later install a newer Boost, configure with:
#     cmake -S . -B build \
#       -DBOOST_ROOT=$HOME/opt/boost-1.87.0 \
#       -DBoost_NO_SYSTEM_PATHS=ON \
#       -DCMAKE_PREFIX_PATH=$HOME/opt/boost-1.87.0

