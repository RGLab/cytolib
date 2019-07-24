
file(GLOB SOURCES "*.cpp")

include_directories(../inst/include)
add_library(cytolib ${SOURCES})
add_dependencies( cytolib h5libs pblibs)

#foreach(dir ${LIBS})
#  message(STATUS "${dir}")
#endforeach()
#static lib can't/doesn't need link
#target_link_libraries(cytolib ${LIBS})

install(TARGETS cytolib DESTINATION lib/cytolib)    