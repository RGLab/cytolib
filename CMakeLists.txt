cmake_minimum_required (VERSION 3.10)
project (cytolib)
# required for unit testing using CMake's ctest command
#    include(CTest)
#    enable_testing()
 
 # required modules for our task
    include(CheckIncludeFile)
    include(CheckIncludeFileCXX)
    include(CheckIncludeFiles)
    include(ExternalProject)
    
    find_package(LAPACK REQUIRED)
    set(LIBS ${LAPACK_LIBRARIES})
    include_directories(${LAPACK_INCLUDE_DIRS})
    
    find_package(BLAS REQUIRED)
    set(LIBS ${LIBS} ${BLAS_LIBRARIES})
    include_directories(${BLAS_INCLUDE_DIRS})

    find_package(Boost REQUIRED)
    include_directories(${Boost_INCLUDE_DIRS})
 
    find_package(ZLIB REQUIRED)
    set(LIBS ${LIBS} ${ZLIB_LIBRARIES})
    include_directories(${ZLIB_INCLUDE_DIRS})
    # TO BE CREATED
#    include(thirdparty/tbb.cmake)
    if ( ${CMAKE_SYSTEM_NAME} STREQUAL "Linux" )
	  	set(LIBS ${LIBS} libdl.so)
  	endif()
#	set(CYTOLIB_INCLUDE_DIR "" CACHE FILEPATH "Path to cytolib includes")
#
#	if(NOT EXISTS ${CYTOLIB_INCLUDE_DIR}/cytolib/global.hpp)
#	  message(SEND_ERROR "Can't find cytolib/global.hpp in ${CYTOLIB_INCLUDE_DIR}")
#	endif()
	
#	include_directories(${CYTOLIB_INCLUDE_DIR})
	  
  # we add the sub-directories that we want CMake to scan
    add_subdirectory(thirdparty)

	include_directories(${THIRDPARTY_INCLUDE_DIR})
#	message(${THIRDPARTY_INCLUDE_DIR})	
	set(CMAKE_CXX_STANDARD 14)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)
	set(LIBS ${LIBS} stdc++fs)

    add_subdirectory(src)
    add_subdirectory(inst)#install header
#    add_subdirectory(test)