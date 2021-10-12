set(EXAMPLES_DIR ${CMAKE_CURRENT_SOURCE_DIR})

include_directories(../tools ../src)

function(example name)
    add_executable(${DIRECTORY_LIB}_${name} ${name}.cpp)
    target_link_libraries(${DIRECTORY_LIB}_${name} PUBLIC world world_tool)
    file(RELATIVE_PATH dir ${EXAMPLES_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
    set_property(GLOBAL APPEND PROPERTY WORLD::EXAMPLES ${dir}/${name})
endfunction(example)

add_subdirectory(analysis_synthesis)
add_subdirectory(codec_test)
add_subdirectory(parameter_io)