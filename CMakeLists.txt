# CMake build system design considerations:
# - Include directories:
#   + Do not define include directories globally using the include_directories
#     command but rather at the target level using the
#     target_include_directories command. That way, it is easier to guarantee
#     that targets are built using the proper list of include directories.
#   + Use the PUBLIC and PRIVATE keywords to specify the scope of include
#     directories. That way, a target linking to a library (using the
#     target_link_librairies command) inherits from the library PUBLIC include
#     directories and not from the PRIVATE ones.

cmake_minimum_required(VERSION 3.0)

# Build Options
option(WORLD_BUILD_TESTS "Set to ON to build tests" OFF)
option(WORLD_BUILD_EXAMPLES "Set to ON to build examples" ON)

# use git version as library version
find_package(Git QUIET)
if (Git_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
    execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-list --count HEAD
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            OUTPUT_VARIABLE git_version
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )
else ()
    set(git_version 0)
endif ()

project(WORLD LANGUAGES CXX VERSION 1.0.${git_version})

add_library(world STATIC
        src/world/cheaptrick.h
        src/world/codec.h
        src/world/common.h
        src/world/constantnumbers.h
        src/world/d4c.h
        src/world/dio.h
        src/world/fft.h
        src/world/harvest.h
        src/world/macrodefinitions.h
        src/world/matlabfunctions.h
        src/world/stonemask.h
        src/world/synthesis.h
        src/world/synthesisrealtime.h
        src/cheaptrick.cpp
        src/codec.cpp
        src/common.cpp
        src/d4c.cpp
        src/dio.cpp
        src/fft.cpp
        src/harvest.cpp
        src/matlabfunctions.cpp
        src/stonemask.cpp
        src/synthesis.cpp
        src/synthesisrealtime.cpp
        )

add_library(world_tool STATIC
        tools/audioio.h
        tools/parameterio.h
        tools/audioio.cpp
        tools/parameterio.cpp
        )

add_library(world::core ALIAS world)
add_library(world::tool ALIAS world_tool)

target_link_libraries(world_tool PUBLIC world)

foreach (lib world world_tool)
    target_include_directories(${lib} PUBLIC $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> PRIVATE src)
    set_target_properties(${lib}
            PROPERTIES
            ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
            )
endforeach ()

include(GNUInstallDirs)
install(TARGETS world world_tool
        EXPORT world-export
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
install(DIRECTORY src/world DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

include(CMakePackageConfigHelpers)
configure_package_config_file(world-config.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/world-config.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/world/cmake
        PATH_VARS CMAKE_INSTALL_INCLUDEDIR
        )
write_basic_package_version_file(
        world-config-version.cmake
        VERSION ${PACKAGE_VERSION}
        COMPATIBILITY AnyNewerVersion
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/world-config.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/world
        )
install(EXPORT world-export
        FILE world-config-version.cmake
        NAMESPACE world::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/world
        )
export(TARGETS world world_tool NAMESPACE world::
        FILE ${PROJECT_BINARY_DIR}/world-targets.cmake)

if (WORLD_BUILD_TESTS)
    add_executable(tests test/test.cpp)
    target_link_libraries(tests
            PRIVATE world
            PRIVATE world_tool
            )
    target_include_directories(tests
            PRIVATE src
            PRIVATE tools
            )
endif ()

if(WORLD_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()