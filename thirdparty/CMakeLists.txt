    add_subdirectory(armadillo)
    
    add_subdirectory(libhdf5)
    
	add_subdirectory(protobuf)
 
	#cp to parent
	set(THIRDPARTY_INCLUDE_DIR ${INCLUDE_DIR} PARENT_SCOPE) 	
	set(LIBS ${H5_CXX_LIBS} ${H5_C_LIBS} ${H5_SZ_LIBS} ${PB_LIBS} ${LIBS} PARENT_SCOPE)
	   