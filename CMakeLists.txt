cmake_minimum_required(VERSION 3.21)
project(episearch LANGUAGES CXX)

# ---------------------------------------------------------------------------
# Standard. C++20 is the target for this exercise: the whole point is to bring
# this code up to something a modern shop would recognise. Do not lower it.
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)   # for clang-tidy / clangd

option(EPISEARCH_WERROR "Treat warnings as errors" OFF)
option(EPISEARCH_TESTS  "Build the test suite"     ON)
option(EPISEARCH_OPENMP "Link OpenMP if available" ON)

# ---------------------------------------------------------------------------
# Warning set. Start with this on; get to zero; then flip EPISEARCH_WERROR on
# and keep it on for the rest of the project.
# ---------------------------------------------------------------------------
add_library(episearch_warnings INTERFACE)
target_compile_options(episearch_warnings INTERFACE
  -Wall
  -Wextra
  -Wpedantic
  -Wshadow                 # a local hiding a member/outer name
  -Wconversion             # silent narrowing, the classic interview trap
  -Wsign-conversion
  -Wold-style-cast         # (int)x instead of static_cast<int>(x)
  -Wcast-align
  -Wunused
  -Woverloaded-virtual
  -Wnon-virtual-dtor
  -Wnull-dereference
  -Wdouble-promotion
  -Wformat=2
  $<$<BOOL:${EPISEARCH_WERROR}>:-Werror>
)

# ---------------------------------------------------------------------------
# Core library: everything except main(). Tests link against this.
# ---------------------------------------------------------------------------
add_library(episearch_core
  srcf/entropy/entropy.cpp
  srcf/utilities/utilities.cpp
  srcf/searchTools/searchTools.cpp
)
target_include_directories(episearch_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/srcf)
target_link_libraries(episearch_core PRIVATE episearch_warnings)

if(EPISEARCH_OPENMP)
  # Apple clang compiles OpenMP but doesn't ship the runtime, and CMake has no
  # idea where Homebrew put it. Without this, find_package quietly fails, the
  # omp pragmas compile to nothing, and you get a serial program that looks
  # like a parallel one.
  if(APPLE)
    execute_process(COMMAND brew --prefix libomp
                    OUTPUT_VARIABLE HOMEBREW_LIBOMP
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET)
    if(HOMEBREW_LIBOMP AND EXISTS "${HOMEBREW_LIBOMP}/lib/libomp.dylib")
      set(OpenMP_ROOT "${HOMEBREW_LIBOMP}")
      set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${HOMEBREW_LIBOMP}/include")
      set(OpenMP_CXX_LIB_NAMES "omp")
      set(OpenMP_omp_LIBRARY "${HOMEBREW_LIBOMP}/lib/libomp.dylib")
    else()
      message(STATUS "OpenMP: no Homebrew libomp found - try `brew install libomp`")
    endif()
  endif()

  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(episearch_core PUBLIC OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP: found, parallel paths enabled")
  else()
    # Apple clang ships without libomp. `brew install libomp` if you want it.
    # Without it the omp pragmas are simply ignored and everything runs serial,
    # which is a perfectly good place to start.
    message(STATUS "OpenMP: NOT found - building serial")
  endif()
endif()

# ---------------------------------------------------------------------------
# Executable
# ---------------------------------------------------------------------------
add_executable(episearch srcf/main.cpp)
target_link_libraries(episearch PRIVATE episearch_core episearch_warnings)

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
if(EPISEARCH_TESTS)
  include(FetchContent)
  FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.15.2
    GIT_SHALLOW    TRUE
  )
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)

  enable_testing()
  add_subdirectory(tests)
endif()
