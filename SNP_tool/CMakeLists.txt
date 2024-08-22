cmake_minimum_required(VERSION 3.11)

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
include(GNUInstallDirs)

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

### Adding ram
find_package(ram 1.0.0 QUIET)
if (NOT ram_FOUND)
    FetchContent_Declare(
            ram
            GIT_REPOSITORY https://github.com/kocijan/ram
            GIT_TAG kocijan/hpc)

    FetchContent_GetProperties(ram)
    if (NOT ram_POPULATED)
        FetchContent_Populate(ram)
        add_subdirectory(
                ${ram_SOURCE_DIR}
                ${ram_BINARY_DIR}
                EXCLUDE_FROM_ALL)
    endif ()
endif ()

# adding edlib
find_package(edlib 1.2.7 QUIET)
if (NOT edlib_FOUND)
    FetchContent_Declare(
            edlib
            GIT_REPOSITORY https://github.com/martinsos/edlib
            GIT_TAG v1.2.7)

    FetchContent_GetProperties(edlib)
    if (NOT edlib_POPULATED)
        FetchContent_Populate(edlib)
        add_subdirectory(
                ${edlib_SOURCE_DIR}
                ${edlib_BINARY_DIR}
                EXCLUDE_FROM_ALL)
    endif ()
endif ()

set_target_properties(ram PROPERTIES POSITION_INDEPENDENT_CODE ON)

include_directories( ${bioparser_INCLUDE_DIRS}  ${biosoup_INCLUDE_DIRS})

add_executable(snp
        src/snp_main.cpp)

# add_library(overlap_source INTERFACE)
# add_library(ram_overlap_source SHARED src/overlap_sources/ram_overlap_source.cpp)
# add_library(overlap_source_factory STATIC src/overlap_sources/overlap_source_factory.cpp)
# add_library(pafEDL_overlap_source SHARED src/overlap_sources/pafEDL_overlap_source.cpp)
# add_library(pafPA_overlap_source SHARED src/overlap_sources/pafPA_overlap_source.cpp)
# add_library(pafM2_overlap_source SHARED src/overlap_sources/pafM2_overlap_source.cpp)
# add_library(pafKSW2_overlap_source SHARED src/overlap_sources/pafKSW2_overlap_source.cpp
#         src/overlap_sources/pafKSW2_overlap_source.cpp)


# target_include_directories(ram_overlap_source PUBLIC
#         ${PROJECT_SOURCE_DIR}/include
#         ${bioparser_INCLUDE_DIRS}
#         ${biosoup_INCLUDE_DIRS}
#         ${thread_pool_INCLUDE_DIRS}
#         ${ram_INCLUDE_DIRS}
#         ${edlib_INCLUDE_DIRS}
#         )


# target_include_directories(overlap_source INTERFACE
#         ${PROJECT_SOURCE_DIR}/include
#         ${bioparser_INCLUDE_DIRS}
#         ${biosoup_INCLUDE_DIRS})

# target_link_libraries(ram_overlap_source
#         bioparser::bioparser
#         biosoup::biosoup
#         thread_pool::thread_pool
#         ram::ram
#         edlib::edlib
#         )

# target_include_directories(pafEDL_overlap_source PUBLIC
#         ${PROJECT_SOURCE_DIR}/include
#         ${bioparser_INCLUDE_DIRS}
#         ${biosoup_INCLUDE_DIRS}
#         )

# target_link_libraries(pafEDL_overlap_source
#         bioparser::bioparser
#         biosoup::biosoup
#         thread_pool::thread_pool
#         edlib::edlib
# )

# target_link_libraries(pafKSW2_overlap_source
#         bioparser::bioparser
#         biosoup::biosoup
#         thread_pool::thread_pool
# )

# target_include_directories( pafPA_overlap_source PUBLIC
#         ${PROJECT_SOURCE_DIR}/include
#         ${biosoup_INCLUDE_DIRS}
# )
# target_link_libraries(pafPA_overlap_source
#         biosoup::biosoup
#         bioparser::bioparser
# )

# target_include_directories( pafM2_overlap_source PUBLIC
#         ${PROJECT_SOURCE_DIR}/include
#         ${biosoup_INCLUDE_DIRS}
# )

# target_link_libraries(pafM2_overlap_source
#         biosoup::biosoup
#         bioparser::bioparser
# )

target_include_directories(snp PUBLIC
        ${PROJECT_SOURCE_DIR}/include
        ${bioparser_INCLUDE_DIRS}
        ${biosoup_INCLUDE_DIRS}
        ${thread_pool_INCLUDE_DIRS}
        ${ram_INCLUDE_DIRS}
        ${edlib_INCLUDE_DIRS}
)

target_link_libraries(snp
        bioparser::bioparser
        biosoup::biosoup
        thread_pool::thread_pool
        ram::ram
        edlib::edlib
        ${CMAKE_DL_LIBS}
        )
