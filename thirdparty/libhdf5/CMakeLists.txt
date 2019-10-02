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
	install (FILES ${HEADERS} DESTINATION include/cytolib)
	
	#bundle into single file
	set(h5_lib_file libhdf5_bundled.a)
	set(MRISCRIPT h5-bundle.mri)
	configure_file(${MRISCRIPT}.in ${MRISCRIPT})
	
	#have to wrap ar cmd into mri script since add_custom_command always prepend space for multi-line cmd
	#which fails makefile that doesn't like space identation
	add_custom_command(OUTPUT ${h5_lib_file}
						COMMAND ar -M < ${MRISCRIPT}
						COMMAND	ranlib ${h5_lib_file} 
						COMMENT "Building the bundled h5 libs")
	add_custom_target(h5libbundle DEPENDS ${h5_lib_file} h5libs)
	set_source_files_properties( ${h5_lib_file} PROPERTIES GENERATED TRUE )
	
	
    
	
    
	set (H5_LIBS_BUNDLE ${CMAKE_CURRENT_BINARY_DIR}/${h5_lib_file} PARENT_SCOPE)
	