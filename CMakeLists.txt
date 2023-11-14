
# =============================================================================

# Minimimum CMake version.
cmake_minimum_required(VERSION 3.14)


# Project name and languages used.
set(PROJECT_NAME caps_sa)
project(${PROJECT_NAME}
        VERSION 0.1.0
        LANGUAGES CXX
)

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Language standards and hard requirements for such.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


# Warning flags.
set(WARNING_FLAGS -Wall -Wextra)

# Extra optimization flags not associated to the `-O` levels.
set(OPTIMIZE_FLAGS -Ofast -mavx2 -march=native)


# Determine the thread library of the system.
include(FindThreads)
if(NOT Threads_FOUND)
    message(FATAL_ERROR "A supported threads library is required. Aborting.")
endif()

# External projects installer.
include(ExternalProject)

set(EXT_INCLUDE ${CMAKE_SOURCE_DIR}/external/include)
file(MAKE_DIRECTORY ${EXT_INCLUDE})


# Prepare the `parlay` library—a parallel algorithms programming toolkit.
message("Build system will fetch and install parlaylib")
ExternalProject_Add(prj_parlaylib
    DOWNLOAD_DIR        ${CMAKE_SOURCE_DIR}/external
    DOWNLOAD_COMMAND    git clone https://github.com/cmuparlay/parlaylib.git
    SOURCE_DIR          ${CMAKE_SOURCE_DIR}/external/parlaylib
    BUILD_IN_SOURCE     TRUE
    INSTALL_DIR         ${CMAKE_SOURCE_DIR}/external
    CONFIGURE_COMMAND   ""
    BUILD_COMMAND       ""
    INSTALL_COMMAND     cp -rf include/parlay ${EXT_INCLUDE}
)


# Default build type.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build Type" FORCE)
endif()


# Add `src` to the build.
add_subdirectory(src)
