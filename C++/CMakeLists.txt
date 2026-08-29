cmake_minimum_required(VERSION 3.13)  # CMake version check
project(straw LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 14)            # Enable c++14 standard
find_package(CURL REQUIRED)
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(ZSTD REQUIRED IMPORTED_TARGET libzstd)

add_library(straw_core STATIC straw.cpp straw_v10.cpp hbs.cpp)
set_target_properties(straw_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(straw_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/straw>)
target_link_libraries(straw_core PUBLIC
  CURL::libcurl ZLIB::ZLIB Threads::Threads PkgConfig::ZSTD)

add_executable(straw_cli main.cpp subsample.cpp compare.cpp)
set_target_properties(straw_cli PROPERTIES OUTPUT_NAME straw)
target_link_libraries(straw_cli PRIVATE straw_core)

include(CTest)

option(STRAW_BUILD_C_API "Build the stable libstraw C API" ON)
if(STRAW_BUILD_C_API)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../C ${CMAKE_CURRENT_BINARY_DIR}/c_api)
endif()

if(BUILD_TESTING)
  find_package(Python3 COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND)
    add_executable(v10_api_probe tests/v10_api_probe.cpp)
    target_include_directories(v10_api_probe PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(v10_api_probe PRIVATE straw_core)
    add_test(NAME v10_reader_conformance
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_v10.py
              $<TARGET_FILE:straw_cli> $<TARGET_FILE:v10_api_probe>)
    add_test(NAME subsample_short
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_subsample.py
              $<TARGET_FILE:straw_cli>)
    add_test(NAME contact_orientation
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_orientation.py
              $<TARGET_FILE:straw_cli> ${CMAKE_CURRENT_SOURCE_DIR}/../R/inst/extdata/test.hic)
  endif()
endif()

option(STRAW_TEST_HTTP "Run tests requiring a loopback HTTP server" OFF)
if(BUILD_TESTING AND Python3_Interpreter_FOUND AND STRAW_TEST_HTTP)
  add_test(NAME v10_http_ranges
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_v10.py
            $<TARGET_FILE:straw_cli> $<TARGET_FILE:v10_api_probe> --http)
endif()
