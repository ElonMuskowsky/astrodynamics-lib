cmake_minimum_required(VERSION 3.16)

project(astrodynamics-lib LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(modules/core)
add_subdirectory(modules/orbital)

enable_testing()
add_subdirectory(tests)

add_executable(astrodynamics-lib main.cpp)

target_link_libraries(astrodynamics-lib PRIVATE core orbital)

include(GNUInstallDirs)
install(TARGETS astrodynamics-lib
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
