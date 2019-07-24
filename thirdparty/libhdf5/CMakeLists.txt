	set(h5_lib_file h5_build/lib/libhdf5_cpp.a)
	file(COPY hdf5small_cxx_1.10.3.tar.gz DESTINATION .)
	file(COPY install.sh DESTINATION .)
	add_custom_command(OUTPUT ${h5_lib_file}
						COMMAND bash "install.sh"
						COMMENT "Building the bundled h5 libs")
	add_custom_target(h5libs DEPENDS ${h5_lib_file})
	set_source_files_properties( ${h5_lib_file} PROPERTIES GENERATED TRUE ) 
	set(h5_include ${CMAKE_CURRENT_BINARY_DIR}/h5_build/include)
	set(INCLUDE_DIR ${INCLUDE_DIR} ${h5_include} PARENT_SCOPE)
	set ( H5_SZ_LIBS ${CMAKE_CURRENT_BINARY_DIR}/h5_build/lib/libsz.a )
	set ( H5_C_LIBS ${CMAKE_CURRENT_BINARY_DIR}/h5_build/lib/libhdf5.a )
	set ( H5_CXX_LIBS ${CMAKE_CURRENT_BINARY_DIR}/h5_build/lib/libhdf5_cpp.a )
	
	#install pb 
	file(GLOB HEADERS "${h5_include}/*")
	install (FILES ${HEADERS} DESTINATION include/cytolib/h5)
	install(FILES ${H5_CXX_LIBS} ${H5_C_LIBS} ${H5_SZ_LIBS} DESTINATION lib/cytolib)  