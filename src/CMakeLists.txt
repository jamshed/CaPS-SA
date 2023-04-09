
# =============================================================================

# Core source files.
set(PROJECT_SRC
        Suffix_Array.cpp
)


# Core library.
add_library(core STATIC ${PROJECT_SRC})


# Include paths for the compiler.
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
target_include_directories(core PUBLIC ${INCLUDE_DIR} ${EXT_INCLUDE})


# Compiler flags.
target_compile_options(core PRIVATE ${WARNING_FLAGS} ${OPTIMIZE_FLAGS})


# Executable for the driver program.
add_executable(${PROJECT_NAME} main.cpp)


# Add dependency of the core to parlaylib.
add_dependencies(core prj_parlaylib)

# Link the core lib to the threads lib.
target_link_libraries(core PRIVATE Threads::Threads)

# Link the executable to the core.
target_link_libraries(${PROJECT_NAME} PRIVATE core)


# Compiler flags for the target.
target_compile_options(${PROJECT_NAME} PRIVATE ${WARNING_FLAGS} ${OPTIMIZE_FLAGS})


# Installation directory.
install(TARGETS ${PROJECT_NAME}
        DESTINATION bin
)

