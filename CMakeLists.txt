cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_COMPILER /usr/local/bin/g++-7)

project(CoMD_RAJA LANGUAGES CXX)

find_package(OpenMP REQUIRED)
if(NOT OPENMP_FOUND)
    message(ERROR "Oh No! No OpenMP Found! :(")
endif()

set(CMAKE_CXX_STANDARD 11)
add_executable(CoMD "")
target_compile_options(CoMD PUBLIC "${OpenMP_CXX_FLAGS}")
target_include_directories(CoMD PUBLIC include)
target_link_libraries(CoMD PUBLIC "${OpenMP_CXX_FLAGS}")

target_sources(CoMD PUBLIC
        src/Atoms.cpp
        src/CoMD.cpp
        src/Command.cpp
        src/Domain.cpp
        src/EAMPotential.cpp
        src/HaloExchange.cpp
        src/LinkCell.cpp
        src/LJPotential.cpp
        src/Parallel.cpp
        src/Random.cpp
        src/Timers.cpp
        src/YAML.cpp)
