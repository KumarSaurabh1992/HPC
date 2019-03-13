cmake_minimum_required(VERSION 3.8)
project(mpi)

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

add_executable(mpi main.cpp common.h common.cpp)
target_link_libraries(mpi ${MPI_LIBRARIES} m)

if(MPI_COMPILE_FLAGS)
    set_target_properties(mpi PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(mpi PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
