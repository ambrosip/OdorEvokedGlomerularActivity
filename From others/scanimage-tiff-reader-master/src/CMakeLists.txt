set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(WIN32)
    file(GLOB PLATFORM_SRCS win32/*.cpp win32/*.h)
else()
    file(GLOB PLATFORM_SRCS linux/*.cpp linux/*.h)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
endif()

find_package(Threads REQUIRED)

# Only the static library is build when using
# the ScanImageTiffReader as a submodule.

set(tgt ScanImageTiffReaderAPIStatic)
SET(SRCS tiff.reader.h
         tiff.reader.api.h
         tiff.reader.api.cpp
         nd.h
    )
add_library(${tgt} STATIC ${SRCS} ${PLATFORM_SRCS})
target_link_libraries(${tgt} Threads::Threads)
target_include_directories(${tgt} PUBLIC ${CMAKE_CURRENT_LIST_DIR})
install(TARGETS ${tgt} DESTINATION lib)
install(FILES tiff.reader.api.h nd.h DESTINATION include)

function(target_add_scanimage_tiff_reader tgt)
    target_link_libraries(${tgt} ScanImageTiffReaderAPIStatic) # should get the include directories
endfunction()

# When the main project is the ScanImageTiffReader
# everything gets built.

if(ScanImageTiffReader STREQUAL ${PROJECT_NAME})
    set(tgt ScanImageTiffReader)
    SET(SRCS main.cpp
             tiff.reader.h
             nd.h
             )
    add_executable(${tgt} ${SRCS} ${PLATFORM_SRCS})
    target_link_libraries(${tgt} Threads::Threads)
    target_include_directories(${tgt} PRIVATE .)
    install(TARGETS ${tgt} DESTINATION bin)

    set(tgt ScanImageTiffReaderAPI)
    SET(SRCS tiff.reader.h
             tiff.reader.api.h
             tiff.reader.api.cpp
             nd.h
        )
    add_library(${tgt} MODULE ${SRCS} ${PLATFORM_SRCS})
    target_link_libraries(${tgt} Threads::Threads)
    target_include_directories(${tgt} PUBLIC .)
    install(TARGETS ${tgt} DESTINATION lib)
ENDIF() # ENDIF ScanImageTiffReader is the current project

## NOTES
#
# 1. In CMake, the install command can only install targets from the current
#    directory (see: https://cmake.org/Bug/view.php?id=14444).
#    In order to make sure the shared library gets installed to the right
#    spots for python and julia, we use little included scripts that end up
#    running from the right directory.
#    At first I had this set up to use add_subdirectory(), but that approach
#    failed due to the above limitation.
#    For matlab, this isn't a problem because the library is statically linked
#    into the mex files.
