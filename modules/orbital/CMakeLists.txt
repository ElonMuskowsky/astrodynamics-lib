add_library(orbital INTERFACE)

target_sources(orbital INTERFACE
    conversions.hpp
    conversions/kepl_cart.hpp
)

target_link_libraries(orbital INTERFACE core)
