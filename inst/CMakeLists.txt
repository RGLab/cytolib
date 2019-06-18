cmake_minimum_required (VERSION 2.6)
project (cytolib)
#target_include_directories(cytolib include/)
file(GLOB HEADERS "include/cytolib/*")

install (FILES ${HEADERS} DESTINATION include/cytolib)

