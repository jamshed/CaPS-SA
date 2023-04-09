
# =============================================================================

# Minimimum CMake version.
cmake_minimum_required(VERSION 3.14)


# Project name and languages used.
set(PROJECT_NAME themis)
project(${PROJECT_NAME}
        VERSION 0.1.0
        LANGUAGES CXX
)


# Language standards and hard requirements for such.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


# Warning flags.
set(WARNING_FLAGS -Wall -Wextra)

# Extra optimization flags not associated to the `-O` levels.
set(OPTIMIZE_FLAGS -funroll-loops)


# Determine the thread library of the system.
include(FindThreads)
if(NOT Threads_FOUND)
    message(FATAL_ERROR "A supported threads library is required. Aborting.")
endif()


# Default build type.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build Type" FORCE)
endif()


# Add `src` to the build.
add_subdirectory(src)
