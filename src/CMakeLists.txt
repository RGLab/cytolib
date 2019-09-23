
file(GLOB SOURCES "*.cpp")

include_directories(../inst/include)
add_library(cytolib SHARED ${SOURCES} GatingSet.pb.cc)#STATIC doesn't work with target_link_libraries
add_dependencies( cytolib h5libs pblibs)
target_link_libraries(cytolib ${LIBS})
set_target_properties(cytolib PROPERTIES VERSION ${CYTOLIB_VERSION}
                                     SOVERSION ${CYTOLIB_VERSION_MAJOR} )
#foreach(dir ${LIBS})
#  message(STATUS "${dir}")
#endforeach()
install(TARGETS cytolib DESTINATION lib)    