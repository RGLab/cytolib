#build cytolib
file(GLOB SOURCES "*.cpp")
include_directories(../inst/include)
add_library(cytolib STATIC ${SOURCES} GatingSet.pb.cc)
add_dependencies( cytolib h5libs pblibs)


install(TARGETS cytolib DESTINATION cytolib/lib)

    