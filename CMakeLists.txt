#
# Copyright (C) 2012-2025 tim cotter. All rights reserved.
#

# search for "you need to"

# you need to set these:
set(THIS_PROJECT_NAME technomancy)
set(THIS_CMAKE_DIR "${CMAKE_SOURCE_DIR}/cmake")

# cmake boilerplate
# unfortunately, we can't put this in preamble.
# see preamble.txt for explanations.
set(THIS_CMAKE_VERSION 3.22.1)
cmake_minimum_required(VERSION ${THIS_CMAKE_VERSION})
message("Creating project ${THIS_PROJECT_NAME}...")
project(${THIS_PROJECT_NAME})

# more cmake boilerplate
include("${THIS_CMAKE_DIR}/scripts/preamble.txt")

# you need to set these:
set(THIS_TARGET_NAME guthrie)

# log it
message("-- Adding executable ${THIS_TARGET_NAME}...")

# gather the source files.
file(GLOB_RECURSE THIS_SOURCE
    ${CMAKE_CURRENT_SOURCE_DIR}/agm/*.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc
)

# gather the header files.
file(GLOB_RECURSE THIS_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/agm/*.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h
)

# add the executable with source and includes in separate groups
add_executable(${THIS_TARGET_NAME} ${THIS_SOURCE} ${THIS_HEADERS})

# add _d to the debug target name
set_target_properties(${THIS_TARGET_NAME} PROPERTIES DEBUG_POSTFIX _d)

# put the binary in bin.
set_target_properties(${THIS_TARGET_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
)

# define the target in the source
target_compile_definitions(${THIS_TARGET_NAME} PRIVATE AGM_TARGET=${THIS_TARGET_NAME})

# add the include directories
set(INCS
    ${CMAKE_SOURCE_DIR}/agm/inc
)
include_directories(${INCS})

# add the libraries
set(LIBS
    curl
)
target_link_libraries(${THIS_TARGET_NAME} ${LIBS})
