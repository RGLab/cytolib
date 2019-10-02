#build cytolib
file(GLOB SOURCES "*.cpp")
include_directories(../inst/include)
add_library(cytolib STATIC ${SOURCES} GatingSet.pb.cc)
add_dependencies( cytolib h5libbundle pblibs)

#bundle all dependencies into single archive
set(CYTOLIB libcytolib.a)
set(cytolib_bundle_file libcytolib-bundled.a)
set(MRISCRIPT cytolib-bundle.mri)
configure_file(${MRISCRIPT}.in ${MRISCRIPT})
add_custom_command(OUTPUT ${cytolib_bundle_file}
						COMMAND ar -M < ${MRISCRIPT}
						COMMAND ranlib ${cytolib_bundle_file}
						COMMENT "Building the bundled cytolib")
add_custom_target(cytolib_bundle ALL DEPENDS ${cytolib_bundle_file} cytolib)
set_source_files_properties( ${cytolib_bundle_file} PROPERTIES GENERATED TRUE ) 

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${cytolib_bundle_file} DESTINATION lib)

    