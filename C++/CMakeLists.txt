cmake_minimum_required(VERSION 3.13)  # CMake version check
project(strawC)               # Create project "simple_example"
set(CMAKE_CXX_STANDARD 14)            # Enable c++14 standard
find_package(CURL REQUIRED)
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(ZSTD REQUIRED IMPORTED_TARGET libzstd)

# Add main.cpp file of project root directory as source file
set(SOURCE_FILES main.cpp straw.cpp straw_v10.cpp subsample.cpp hbs.cpp compare.cpp)
add_executable(straw ${SOURCE_FILES})

target_link_libraries(straw CURL::libcurl ZLIB::ZLIB Threads::Threads PkgConfig::ZSTD)

include(CTest)
if(BUILD_TESTING)
  find_package(Python3 COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND)
    add_executable(v10_api_probe tests/v10_api_probe.cpp straw.cpp straw_v10.cpp)
    target_include_directories(v10_api_probe PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(v10_api_probe PRIVATE CURL::libcurl ZLIB::ZLIB Threads::Threads PkgConfig::ZSTD)
    add_test(NAME v10_reader_conformance
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_v10.py
              $<TARGET_FILE:straw> $<TARGET_FILE:v10_api_probe>)
    add_test(NAME subsample_short
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_subsample.py
              $<TARGET_FILE:straw>)
  endif()
endif()

option(STRAW_TEST_HTTP "Run tests requiring a loopback HTTP server" OFF)
if(BUILD_TESTING AND Python3_Interpreter_FOUND AND STRAW_TEST_HTTP)
  add_test(NAME v10_http_ranges
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_v10.py
            $<TARGET_FILE:straw> $<TARGET_FILE:v10_api_probe> --http)
endif()
