add_library(core INTERFACE)

target_sources(core INTERFACE
    math.hpp
    constants.hpp
    state.hpp
)

target_include_directories(core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/..)
