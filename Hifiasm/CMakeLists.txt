cmake_minimum_required(VERSION 3.16.3)

project(hifiasm
        VERSION 1.0
        DESCRIPTION ""
        LANGUAGES CXX
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

include(FetchContent)

## Find bioparser package
find_package(bioparser 3.0.13 QUIET)
if (NOT bioparser_FOUND)
    FetchContent_Declare(
            bioparser
            GIT_REPOSITORY https://github.com/rvaser/bioparser
            GIT_TAG 3.0.13)
    FetchContent_GetProperties(bioparser)
    if (NOT bioparser_POPULATED)
        FetchContent_Populate(bioparser)
        add_subdirectory(
                ${bioparser_SOURCE_DIR}
                ${bioparser_BINARY_DIR}
                EXCLUDE_FROM_ALL)
    endif ()
endif ()

## Find biosoup package
find_package(biosoup 0.11.0 QUIET)
if (NOT bioparser_FOUND)
    FetchContent_Declare(
            biosoup
            GIT_REPOSITORY https://github.com/rvaser/biosoup
            GIT_TAG 0.11.0)
    FetchContent_GetProperties(biosoup)
    if (NOT biosoup_POPULATED)
        FetchContent_Populate(biosoup)
        add_subdirectory(
                ${biosoup_SOURCE_DIR}
                ${biosoup_BINARY_DIR}
                EXCLUDE_FROM_ALL)
    endif ()
endif ()

## Adding tests

find_package(GTest 1.11.0 QUIET)
if (NOT GTest_FOUND)
    FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest
            GIT_TAG release-1.11.0)

    FetchContent_GetProperties(googletest)
    if (NOT googletest_POPULATED)
        FetchContent_Populate(googletest)
        add_subdirectory(
                ${googletest_SOURCE_DIR}
                ${googletest_BINARY_DIR}
                EXCLUDE_FROM_ALL)
        add_library(GTest::Main ALIAS gtest_main)
    endif ()
endif ()

### Adding thread pool
find_package(threadPool 4.0.0 QUIET)
if (NOT thread_pool_FOUND)
    FetchContent_Declare(
            thread_pool
            GIT_REPOSITORY https://github.com/rvaser/thread_pool
            GIT_TAG 4.0.0)
    FetchContent_GetProperties(thread_pool)
    if (NOT thread_pool_POPULATED)
        FetchContent_Populate(thread_pool)
        add_subdirectory(
                ${thread_pool_SOURCE_DIR}
                ${thread_pool_BINARY_DIR}
                EXCLUDE_FROM_ALL)
    endif ()
endif ()

add_executable(hifiasm_exe
src/hifiasm.cpp)

target_link_libraries(hifiasm_exe
        bioparser::bioparser
        biosoup::biosoup
        thread_pool::thread_pool
)