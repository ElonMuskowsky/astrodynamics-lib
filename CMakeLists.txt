cmake_minimum_required(VERSION 3.16)

project(astrodynamics-lib LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(astrodynamics-lib main.cpp)

include(GNUInstallDirs)
install(TARGETS astrodynamics-lib
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
