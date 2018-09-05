# Setting up external library for FASTAHACK
SET(FASTAHACK_PROJECT fastahack_project CACHE INTERNAL "fastahack project name")
SET(FASTAHACK_DIR ${CMAKE_BINARY_DIR}/externals/fastahack CACHE INTERNAL "fastahack project directory")
ExternalProject_Add(${FASTAHACK_PROJECT}
	GIT_REPOSITORY https://github.com/dillonl/fastahack.git
	GIT_TAG master
	INSTALL_COMMAND ""
    UPDATE_COMMAND ""
    PREFIX ${FASTAHACK_DIR}
    CMAKE_CACHE_ARGS
        -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}

)

ExternalProject_Get_Property(${FASTAHACK_PROJECT} INSTALL_DIR)
ExternalProject_Get_Property(${FASTAHACK_PROJECT} SOURCE_DIR)
ExternalProject_Get_Property(${FASTAHACK_PROJECT} BINARY_DIR)

SET(FASTAHACK_LIB ${BINARY_DIR}/src/libfastahack_src_lib.a CACHE INTERNAL "FASTAHACK Library")
SET(FASTAHACK_INCLUDE ${SOURCE_DIR}/src CACHE INTERNAL "FASTAHACK Include")