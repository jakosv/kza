cmake_minimum_required(VERSION 3.10)

project(KZ)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -fPIC")

add_library(${PROJECT_NAME} SHARED kz.cpp kza.cpp)

if (kza_timer)
    add_definitions(-DTIMER)
endif()
if (prefix_sum)
    add_definitions(-DPREFIX_SUM)
endif()
if (NOT no_python)
    set(PYBIND11_FINDPYTHON ON)
    set(PYBIND11_PYTHON_VERSION 3.6)
    find_package(pybind11 REQUIRED)

    pybind11_add_module(libKZ_py pywrap.cpp kz.cpp kza.cpp)

    target_include_directories(${PROJECT_NAME} PRIVATE ${PYBIND11_INCLUDE_DIRS})

    set_target_properties(libKZ_py
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "../build"
    LIBRARY_OUTPUT_DIRECTORY "../build"
    RUNTIME_OUTPUT_DIRECTORY "../build"
    )
endif()
    
target_compile_definitions(${PROJECT_NAME} PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})

set_target_properties(${PROJECT_NAME}
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "../build"
    LIBRARY_OUTPUT_DIRECTORY "../build"
    RUNTIME_OUTPUT_DIRECTORY "../build"
)