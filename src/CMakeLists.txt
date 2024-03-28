cmake_minimum_required(VERSION 3.10)

project(KZ_Lib)
set(CMAKE_CXX_STANDARD 20)
set(PYBIND11_PYTHON_VERSION 3.6)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -fPIC")

find_package(pybind11 REQUIRED)

if (kza_timer)
    add_definitions(-DTIMER)
endif()
if (prefix_sum)
    add_definitions(-DPREFIX_SUM)
endif()

pybind11_add_module(${PROJECT_NAME} pywrap.cpp kz.cpp kza.c timer.c)

target_compile_definitions(${PROJECT_NAME} PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
target_include_directories(${PROJECT_NAME} PRIVATE ${PYBIND11_INCLUDE_DIRS})

set_target_properties( KZ_Lib
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "../tests"
    LIBRARY_OUTPUT_DIRECTORY "../tests"
    RUNTIME_OUTPUT_DIRECTORY "../tests"
)