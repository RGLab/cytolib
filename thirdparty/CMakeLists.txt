    add_subdirectory(armadillo)
    
    add_subdirectory(libhdf5)
    
	#cp to parent
	set(THIRDPARTY_INCLUDE_DIR ${INCLUDE_DIR} PARENT_SCOPE) 	
#	set ( PB_LIBS ${PB_LIBS} PARENT_SCOPE)
	   