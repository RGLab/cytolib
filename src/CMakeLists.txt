#Copyright 2019 Fred Hutchinson Cancer Research Center
#See the included LICENSE file for details on the licence that is granted to the user of this software.
#build cytolib
file(GLOB SOURCES "*.cpp")
include_directories(../inst/include)
add_library(cytolib STATIC ${SOURCES} GatingSet.pb.cc)
add_dependencies( cytolib h5libs pblibs)

find_package(TileDB)

if(TILEDB_FOUND)
add_definitions(-DHAVE_TILEDB)
#include_directories(${TILEDB_LIB}/include)
target_link_libraries(cytolib TileDB::tiledb_shared)
endif()

install(TARGETS cytolib DESTINATION cytolib/lib)

    
