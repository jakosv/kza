cmake_minimum_required(VERSION 3.10)

project(KZ)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-Wall -fPIC")

if (KZA_TIMER)
    add_definitions(-DKZA_TIMER)
endif()
if (PREFIX_SUM)
    add_definitions(-DPREFIX_SUM)
endif()
set(PYBIND11_FINDPYTHON ON)
set(PYBIND11_PYTHON_VERSION 3.6)
find_package(pybind11 REQUIRED)

pybind11_add_module(libKZ_py pywrap.cpp)

set_target_properties(libKZ_py
PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "."
    LIBRARY_OUTPUT_DIRECTORY "."
    RUNTIME_OUTPUT_DIRECTORY "."
)
