cmake_minimum_required(VERSION 2.6)

project(soap)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

include(conan)
conan_cmake_run(CONANFILE conanfile.py
                BASIC_SETUP CMAKE_TARGETS
                BUILD missing)


add_executable(soap soap.cpp)
target_link_libraries(soap ${CMAKE_DL_LIBS} CONAN_PKG::boost CONAN_PKG::VTK-maxi CONAN_PKG::eigen)