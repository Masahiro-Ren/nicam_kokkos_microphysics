cmake_minimum_required(VERSION 3.16)
project(NICAM_DCKERNEL LANGUAGES CXX)

# include architecture settings
include(config/archi_flags.cmake)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Kokkos (assumes it was installed with CMake)
# !!!!! Dont forget Specify -DKokkos_ROOT first !!!!!!!!
find_package(Kokkos REQUIRED)

# Include header files
include_directories(include)

# Collect all source files from src/
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS src/*.cpp)

# Create executable
add_executable(phy.exe ${SOURCES})

# Link Kokkos
target_link_libraries(phy.exe Kokkos::kokkos)

# Optional: show compile flags
message(STATUS "Compiling with CXX Flags: ${CMAKE_CXX_FLAGS}")
