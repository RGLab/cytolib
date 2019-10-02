    add_subdirectory(armadillo)
    
    add_subdirectory(libhdf5)
    
	add_subdirectory(protobuf)
 
	#cp to parent
	set(THIRDPARTY_INCLUDE_DIR ${INCLUDE_DIR} PARENT_SCOPE) 	
	set (H5_LIBS_BUNDLE ${H5_LIBS_BUNDLE} PARENT_SCOPE)
	set ( PB_LIBS ${PB_LIBS} PARENT_SCOPE)
	   