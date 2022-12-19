cmake_minimum_required (VERSION 3.12)

project(LocationRouting)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${CMAKE_CURRENT_SOURCE_DIR}/../../Bapcod/CMake ${CMAKE_CURRENT_SOURCE_DIR}/../../Bapcod/CMake/Packages
		      ${CMAKE_CURRENT_SOURCE_DIR}/CMake ${CMAKE_CURRENT_SOURCE_DIR}/CMake/Packages)

find_package(BcpRcsp)
find_package(NonPublicCuts)

include(BapcodApplicationProject)
include(BapcodApplicationPackage)
include(Documentation)
include(Profile)

set(USER_LIBRARIES ${BCP_RCSP_LIBRARY})
include_directories(${BCP_RCSP_INCLUDE_DIR})

# Execute the build
set_project_version(1 0 0)

# Build Application
try_build_bapcod_application()

# Make Package
try_build_bapcod_application_package()
