
file(GLOB SOURCES "*.cpp")

include_directories(../inst/include)
add_library(cytoliblite ${SOURCES})
add_dependencies( cytoliblite h5libs pblibs)

#foreach(dir ${LIBS})
#  message(STATUS "${dir}")
#endforeach()
#install(TARGETS cytolib DESTINATION lib)    