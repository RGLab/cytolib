
file(GLOB SOURCES "*.cpp")

include_directories(../inst/include)
add_library(cytolib ${SOURCES})
add_dependencies( cytolib h5libs pblibs)

#foreach(dir ${LIBS})
#  message(STATUS "${dir}")
#endforeach()
install(TARGETS cytolib DESTINATION lib)    