cmake_minimum_required(VERSION 2.6)
set(CMAKE_CXX_STANDARD 14)

project(soap)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

include(conan)
conan_cmake_run(CONANFILE conanfile.py
                BASIC_SETUP CMAKE_TARGETS
                BUILD missing)

file(GLOB SOAPFILES ${PROJECT_SOURCE_DIR}/*.h ${PROJECT_SOURCE_DIR}/*.cpp)
add_executable(soap ${SOAPFILES})
target_link_libraries(soap ${CMAKE_DL_LIBS} CONAN_PKG::boost CONAN_PKG::VTK-maxi CONAN_PKG::eigen)
