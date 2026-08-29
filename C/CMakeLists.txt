include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

add_library(straw_c_api SHARED src/straw_c.cpp)
set_target_properties(straw_c_api PROPERTIES
  OUTPUT_NAME straw
  EXPORT_NAME straw
  VERSION 1.0.0
  SOVERSION 1
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN YES)
target_compile_definitions(straw_c_api PRIVATE STRAW_BUILDING_LIBRARY)
target_include_directories(straw_c_api PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(straw_c_api PRIVATE straw_core)
add_library(straw::c_api ALIAS straw_c_api)

configure_package_config_file(
  cmake/strawConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/strawConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/straw)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/strawConfigVersion.cmake
  VERSION 1.0.0
  COMPATIBILITY SameMajorVersion)
configure_file(libstraw.pc.in ${CMAKE_CURRENT_BINARY_DIR}/libstraw.pc @ONLY)

install(TARGETS straw_c_api
  EXPORT strawTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT strawTargets
  FILE strawTargets.cmake
  NAMESPACE straw::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/straw)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/strawConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/strawConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/straw)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libstraw.pc
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

if(BUILD_TESTING)
  add_executable(straw_c_smoke tests/test_straw_c.c)
  set_target_properties(straw_c_smoke PROPERTIES LINKER_LANGUAGE CXX)
  target_link_libraries(straw_c_smoke PRIVATE straw_c_api)
  add_test(NAME c_api_smoke
    COMMAND straw_c_smoke ${CMAKE_CURRENT_SOURCE_DIR}/../R/inst/extdata/test.hic)
endif()
