cmake_minimum_required(VERSION 2.8)
project(World)

set(PROJECT_VER_MAJOR 0)
set(PROJECT_VER_MINOR 3)
set(PROJECT_VER_PATCH 0)
set(PROJECT_VER "${PROJECT_VER_MAJOR}.${PROJECT_VER_MINOR}.${PROJECT_VER_PATCH}")
set(PROJECT_APIVER "${PROJECT_VER_MAJOR}.${PROJECT_VER_MINOR}")

option(BUILD_EXAMPLES "Build examples" ON)

# set the default path for built executables to the "bin" directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# set the default path for built libraries to the "lib" directory
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

# generate header
set(CONFIG_FILE "${PROJECT_BINARY_DIR}/world/config.h")
configure_file("src/world/config.h.in" ${CONFIG_FILE} @ONLY)

configure_file(world.pc.in "${PROJECT_BINARY_DIR}/world.pc" @ONLY)

file(GLOB SOURCES src/*.cpp)
file(GLOB INCLUDES src/world/*.h)

include_directories(
  "src"
  ${PROJECT_BINARY_DIR} # for generated headers
  )

# Core library
add_library(world ${SOURCES} ${INCLUDES} ${CONFIG_FILE})
set_target_properties(world PROPERTIES
  # TODO
  # CXX_VISIBILITY_PRESET hidden
  # VISIBILITY_INLINES_HIDDEN 1
  VERSION ${PROJECT_VER}
  SOVERSION ${PROJECT_APIVER}
)

install(TARGETS world DESTINATION lib RUNTIME DESTINATION bin)
install(DIRECTORY src/world DESTINATION include PATTERN "*.in" EXCLUDE)
install(DIRECTORY "${PROJECT_BINARY_DIR}/world" DESTINATION include)
install(FILES "${PROJECT_BINARY_DIR}/world.pc" DESTINATION lib/pkgconfig/)
