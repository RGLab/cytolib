	file(COPY hdf5small_cxx_1.10.3.tar.gz DESTINATION .)
	file(COPY install.sh DESTINATION .)
	set ( H5_SZ_LIBS ${CMAKE_CURRENT_BINARY_DIR}/h5_build/lib/libsz.a)
	set ( H5_C_LIBS ${CMAKE_CURRENT_BINARY_DIR}/h5_build/lib/libhdf5.a)
	set ( H5_CXX_LIBS ${CMAKE_CURRENT_BINARY_DIR}/h5_build/lib/libhdf5_cpp.a)
	add_custom_command(OUTPUT ${H5_CXX_LIBS}
						COMMAND bash "install.sh"
						COMMENT "Building h5")
	add_custom_target(h5libs DEPENDS ${H5_CXX_LIBS})
	set_source_files_properties( ${H5_CXX_LIBS} PROPERTIES GENERATED TRUE )
	set(h5_include ${CMAKE_CURRENT_BINARY_DIR}/h5_build/include)
    set(INCLUDE_DIR ${INCLUDE_DIR} ${h5_include} PARENT_SCOPE)
	file(GLOB HEADERS "${h5_include}/*")
	install (FILES ${HEADERS} DESTINATION cytolib/include)
	install (FILES ${H5_SZ_LIBS} DESTINATION cytolib/lib)
	install (FILES ${H5_C_LIBS} DESTINATION cytolib/lib)
	install (FILES ${H5_CXX_LIBS} DESTINATION cytolib/lib)
	
		