# Copyright (c) 2017 Vidrio Technologies, All Rights Reserved
# Author: Nathan Clack <nathan@vidriotech.com>
find_program(GIT git)
if(GIT)
    exec_program(${GIT} ${CMAKE_CURRENT_SOURCE_DIR} ARGS "describe --tags --abbrev=0" OUTPUT_VARIABLE GIT_TAG)
    exec_program(${GIT} ${CMAKE_CURRENT_SOURCE_DIR} ARGS "describe --always" OUTPUT_VARIABLE GIT_HASH)
    add_definitions(-DGIT_TAG=${GIT_TAG})
    add_definitions(-DGIT_HASH=${GIT_HASH})
    set(CPACK_PACKAGE_VERSION ${GIT_TAG})
    message("Version ${GIT_TAG} ${GIT_HASH}")
else()
    #add_definitions(-DGIT_TAG="Unknown" -DGIT_HASH=" ")
endif()
