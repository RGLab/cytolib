	set(pb_lib_file pb_build/lib/libprotobuf-lite.a)
	file(COPY protobuf-3.7.0.tar.gz DESTINATION .)
	file(COPY install.sh DESTINATION .)
	add_custom_command(OUTPUT ${pb_lib_file}
						COMMAND bash install.sh
						COMMENT "Building the bundled protobuf libs")
	add_custom_target(pblibs DEPENDS ${pb_lib_file})
	set_source_files_properties( ${pb_lib_file} PROPERTIES GENERATED TRUE ) 
	
   	set(INCLUDE_DIR ${INCLUDE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/pb_build/include PARENT_SCOPE)
	set ( PB_LIBS ${CMAKE_CURRENT_BINARY_DIR}/pb_build/lib/libprotobuf-lite.a PARENT_SCOPE)
