	set(pb_lib_file pb_build/lib/libprotobuf-lite.a)
	file(COPY protobuf-3.10.0.tar.gz DESTINATION .)
	file(COPY install.sh DESTINATION .)
	add_custom_command(OUTPUT ${pb_lib_file}
						COMMAND bash install.sh
						COMMENT "Building the protobuf libs")
	add_custom_target(pblibs DEPENDS ${pb_lib_file})
	set_source_files_properties( ${pb_lib_file} PROPERTIES GENERATED TRUE ) 
	
    set(pb_include ${CMAKE_CURRENT_BINARY_DIR}/pb_build/include)
    set(INCLUDE_DIR ${INCLUDE_DIR} ${pb_include} PARENT_SCOPE)
    
    install (DIRECTORY "${pb_include}/google" DESTINATION cytolib/include)
    install (FILES ${CMAKE_CURRENT_BINARY_DIR}/pb_build/lib/libprotobuf-lite.a DESTINATION cytolib/lib)
