cmake_minimum_required (VERSION 3.10)
project (cytolib)
# required for unit testing using CMake's ctest command
#    include(CTest)
#    enable_testing()
	set(CYTOLIB_VERSION_MAJOR 2)
	set(CYTOLIB_VERSION_MINOR 2)
	set(CYTOLIB_VERSION_PATCH 17)
	set(CYTOLIB_VERSION ${CYTOLIB_VERSION_MAJOR}.${CYTOLIB_VERSION_MINOR}.${CYTOLIB_VERSION_PATCH} )
 # required modules for our task
    include(CheckIncludeFile)
    include(CheckIncludeFileCXX)
    include(CheckIncludeFiles)
    include(ExternalProject)
    
	  
  # we add the sub-directories that we want CMake to scan
    add_subdirectory(thirdparty)

	include_directories(${THIRDPARTY_INCLUDE_DIR})
	set(CMAKE_CXX_STANDARD 14)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)

    add_subdirectory(src)
    add_subdirectory(inst)#install header
#    add_subdirectory(test)