#
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
#
# CMake options for build WORLD
# - STATIC_LIB: For build the static library, default off
#

cmake_minimum_required(VERSION 3.0)
project(WORLD CXX)

# options
option(STATIC_LIB "Build as static library" OFF)

file(GLOB SRCS "src/*.cpp" "src/*.c")
file(GLOB HEADERS "src/world/*.h")


if (STATIC_LIB)
    add_library(${PROJECT_NAME} ${SRCS} ${HEADERS})
else (STATIC_LIB)
    add_library(${PROJECT_NAME} SHARED ${SRCS} ${HEADERS})
endif (STATIC_LIB)

target_include_directories(${PROJECT_NAME} PUBLIC "src")

install(TARGETS ${PROJECT_NAME})
