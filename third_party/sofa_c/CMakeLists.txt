file(GLOB_RECURSE SOFA_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.c")
file(GLOB_RECURSE _sofa_h "${CMAKE_CURRENT_SOURCE_DIR}/sofa.h")

list(GET _sofa_h 0 _sofa_h_path)
get_filename_component(SOFA_INCLUDE_DIR "${_sofa_h_path}" DIRECTORY)

add_library(sofa STATIC ${SOFA_SOURCES})
target_include_directories(sofa PUBLIC "${SOFA_INCLUDE_DIR}")
