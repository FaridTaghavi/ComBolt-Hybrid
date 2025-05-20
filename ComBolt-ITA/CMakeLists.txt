# CMakeLists.txt

cmake_minimum_required(VERSION 3.10)
project(ComBolt-ITA)


set(USE_ROOT TRUE)

if (USE_ROOT)
    message(STATUS "Looking for ROOT")
    find_package(ROOT 6.21 CONFIG REQUIRED)
    message(STATUS "ROOT ${ROOT_VERSION} found at ${ROOT_DIR}.")
    include(${ROOT_USE_FILE})
endif ()

# target_link_libraries(kineticTheory /usr/lib/x86_64-linux-gnu/libpthread.a)



# find_library(PTHREAD_LIBRARY pthread)
# target_link_libraries(kineticTheory ${PTHREAD_LIBRARY})
# target_link_libraries(kineticTheory pthread c)


# message(STATUS "Looking for BOOST library")
# find_package(Boost COMPONENTS filesystem program_options system REQUIRED)
# message(STATUS "BOOST include dir found at ${Boost_INCLUDE_DIR}")
#link_directories(${Boost_LIBRARY_DIR})



message(STATUS "Looking for BOOST library")

# Set custom paths for Boost
# set(BOOST_ROOT "/usr/local")
# set(BOOST_INCLUDEDIR "/usr/local/include")
# set(BOOST_LIBRARYDIR "/usr/local/lib")

# list(APPEND CMAKE_PREFIX_PATH "/usr/local")

message(STATUS "Looking for BOOST library in ${BOOST_ROOT}")

# Now find Boost
find_package(Boost 1.87 REQUIRED COMPONENTS filesystem program_options system math)

message(STATUS "Boost include directories: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost library directories: ${Boost_LIBRARY_DIRS}")
# link_directories(${Boost_LIBRARY_DIR})

# Find HDF5 package
message(STATUS "Looking for HDF5 library")
find_package(HDF5 COMPONENTS C CXX REQUIRED)
# find_package(HDF5 CONFIG REQUIRED)
message(STATUS "HDF5 include directory: ${HDF5_INCLUDE_DIRS}")
message(STATUS "HDF5 libraries: ${HDF5_LIBRARIES}")

# Find GSL package
message(STATUS "Looking for GSL library")
find_package(GSL REQUIRED)
message(STATUS "GSL include dir found at ${GSL_INCLUDE_DIR}")

# Find OpenMP
message(STATUS "Looking for OpenMP library")
find_package(OpenMP REQUIRED)
message(STATUS "OpenMP OpenMP found at ${OpenMP_INCLUDE_DIR}")


include_directories(${CMAKE_SOURCE_DIR}/src SYSTEM ${Boost_INCLUDE_DIR} ${HDF5_INCLUDE_DIRS} ${GSL_INCLUDE_DIR})
add_subdirectory(src)



set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++17" )

set(CMAKE_EXE_LINKER_FLAGS "-L/usr/lib/x86_64-linux-gnu")


# Add the main file
add_executable(ComBolt-ITA ComBolt-ITA.cpp)
target_link_libraries(ComBolt-ITA src ${ROOT_LIBRARIES} ${Boost_LIBRARIES} stdc++fs ${HDF5_LIBRARIES} GSL::gsl GSL::gslcblas OpenMP::OpenMP_CXX)

# target_link_libraries(kineticTheory src ${ROOT_LIBRARIES} Boost::Math ${Boost_LIBRARIES} stdc++fs ${HDF5_LIBRARIES} GSL::gsl GSL::gslcblas OpenMP::OpenMP_CXX)
# target_link_libraries(kineticTheory src ${ROOT_LIBRARIES} Boost::Math ${Boost_LIBRARIES} stdc++fs ${HDF5_LIBRARIES} GSL::gsl GSL::gslcblas OpenMP::OpenMP_CXX )


# target_link_libraries(kineticTheory /usr/lib/x86_64-linux-gnu/libpthread.a /usr/lib/x86_64-linux-gnu/libc.a)

# Add source files
#target_sources(kineticTheory PRIVATE
#    #src/functions.cpp
#    #src/administration.cpp
#    #src/initial_state.cpp
#    #src/rungeKutta_Solver.cpp
#    src/RK4MoC.cpp
#)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -fsanitize=address")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -g" )

#target_link_libraries(kineticTheory PRIVATE ${HDF5_LIBRARIES})

