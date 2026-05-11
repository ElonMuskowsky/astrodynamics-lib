include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP ON
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(tst_kepl_cart tst_kepl_cart.cpp)
target_link_libraries(tst_kepl_cart PRIVATE GTest::gtest_main orbital)

include(GoogleTest)
gtest_discover_tests(tst_kepl_cart)