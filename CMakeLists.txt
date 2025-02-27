
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


if(ASAN OR UBSAN)
    set(CMAKE_C_COMPILER clang)
    set(CMAKE_CXX_COMPILER clang++)
    # set(CMAKE_LINKER_TYPE LLD)
endif()

# Warning flags.
set(WARNING_FLAGS -Wall -Wextra -Wpedantic)
set(SUPPRESS_WARNING_FLAGS -Wno-self-assign-overloaded)

# Bundle debugging information flags that we want to conditionally
# pass on to the compiler.
# Reference: https://gcc.gnu.org/onlinedocs/gcc/Debugging-Options.html.
set(DEBUG_FLAGS -ggdb)

# Bundle address sanitization flags: https://github.com/google/sanitizers/wiki/AddressSanitizer.
set(ASAN_FLAGS -fsanitize=address -fno-omit-frame-pointer)
# set(ASAN_LINK_FLAGS -static-libasan)

# Bundle undefined behavior sanitization flags: https://clang.llvm.org/docs/UndefinedBehaviorSanitizer.html.
set(UBSAN_FLAGS -fsanitize=undefined
                -fsanitize=float-divide-by-zero #-fsanitize=unsigned-integer-overflow
                -fsanitize=implicit-conversion -fsanitize=local-bounds -fsanitize=nullability)
# set(UBSAN_LINK_FLAGS -fuse-ld=lld)

# Extra optimization flags not associated to the `-O` levels.
set(OPTIMIZE_FLAGS -Ofast -mavx2 -funroll-loops -march=native)


# Accumulate the compile- and link-flags.
set(COMPILE_FLAGS "")
set(LINK_FLAGS "")

# Specify the warnings and the extra optimization flags to the compiler for the target library.
list(APPEND COMPILE_FLAGS ${WARNING_FLAGS} ${SUPPRESS_WARNING_FLAGS} ${OPTIMIZE_FLAGS})

if(ASAN)
    list(APPEND COMPILE_FLAGS ${ASAN_FLAGS} ${DEBUG_FLAGS})
    list(APPEND LINK_FLAGS ${ASAN_LINK_FLAGS} ${ASAN_FLAGS})
endif()

if(UBSAN)
    list(APPEND COMPILE_FLAGS ${UBSAN_FLAGS} ${DEBUG_FLAGS})
    list(APPEND LINK_FLAGS ${UBSAN_LINK_FLAGS} ${UBSAN_FLAGS})
endif()


# Add the required preprocessor definitions (`#define`s) to pass on.
add_compile_definitions(CMAKE_EXPORT_COMPILE_COMMANDS)  # For better `clangd` scanning.


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
