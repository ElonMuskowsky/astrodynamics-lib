find_package(Qt6 REQUIRED COMPONENTS Test)

qt_add_executable(tst_kepl_cart tst_kepl_cart.cpp)

target_link_libraries(tst_kepl_cart PRIVATE Qt6::Test orbital)

add_test(NAME tst_kepl_cart COMMAND tst_kepl_cart)
