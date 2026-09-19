add_executable(episearch_tests
  test_entropy.cpp
)

target_link_libraries(episearch_tests
  PRIVATE
    episearch_core
    episearch_warnings
    GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(episearch_tests)
