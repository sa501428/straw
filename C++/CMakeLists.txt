cmake_minimum_required(VERSION 3.13)  # CMake version check
project(strawC)               # Create project "simple_example"
set(CMAKE_CXX_STANDARD 14)            # Enable c++14 standard
find_package(CURL REQUIRED)
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(ZSTD REQUIRED IMPORTED_TARGET libzstd)

# Add main.cpp file of project root directory as source file
set(SOURCE_FILES main.cpp straw.cpp)
add_executable(straw ${SOURCE_FILES})

target_link_libraries(straw CURL::libcurl ZLIB::ZLIB Threads::Threads PkgConfig::ZSTD)
