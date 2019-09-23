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
	if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
		set(FSLIB c++fs)
	elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  	  	set(FSLIB stdc++fs)
	else()
	  	message(SEND_ERROR "unsupported compiler ${CMAKE_CXX_COMPILER_ID}")
	endif()
#  	message("compiler ${CMAKE_CXX_COMPILER_ID}")
	set(LIBS ${LIBS} ${FSLIB})

    add_subdirectory(src)
    add_subdirectory(inst)#install header
#    add_subdirectory(test)