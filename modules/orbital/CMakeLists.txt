add_library(orbital INTERFACE)

target_sources(orbital INTERFACE
    kepl_cart_vallado.hpp
)

target_link_libraries(orbital INTERFACE core)
