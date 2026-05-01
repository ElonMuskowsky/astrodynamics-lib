add_library(frames STATIC frame_conversions.cpp)

target_include_directories(frames PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(frames
    PUBLIC  core
    PRIVATE sofa
)
