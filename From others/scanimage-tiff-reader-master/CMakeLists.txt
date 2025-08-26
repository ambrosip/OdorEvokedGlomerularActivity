cmake_minimum_required (VERSION 3.0)
project(ScanImageTiffReader)

include(cmake/git-versioning.cmake)
include(cmake/static-runtime.cmake)

set_property(GLOBAL PROPERTY USE_FOLDERS TRUE)
set(CMAKE_C_STANDARD   99)
set(CMAKE_CXX_STANDARD 11)

enable_testing()
add_subdirectory(src)
#add_subdirectory(data) # runs some tests over this data
#add_subdirectory(tests)

include(CPack) # enable packaging
