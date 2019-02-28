
SET(SAMTOOLS_INSTALL make prefix=${PROJECT_SOURCE_DIR}/bin/externals/samtools/src/samtools_project install) 

SET(SAMTOOLS_PROJECT samtools_project CACHE INTERNAL "samtools project name")
SET(SAMTOOLS_DIR ${CMAKE_BINARY_DIR}/externals/samtools CACHE INTERNAL "samtools project directory")

ExternalProject_Add(${SAMTOOLS_PROJECT}
        URL  https://github.com/samtools/samtools/archive/1.9.tar.gz
	DEPENDS ${HTSLIB_PROJECT}
        CONFIGURE_COMMAND ""
        BUILD_IN_SOURCE 1
    	BUILD_COMMAND make HTSDIR=${HTSLIB_INCLUDE_DIR}
   	INSTALL_COMMAND "${SAMTOOLS_INSTALL}"
        UPDATE_COMMAND ""
        PREFIX ${SAMTOOLS_DIR}
    	CMAKE_CACHE_ARGS
		-DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
        	-DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}
)

ExternalProject_Get_Property(${SAMTOOLS_PROJECT} INSTALL_DIR)
ExternalProject_Get_Property(${SAMTOOLS_PROJECT} SOURCE_DIR)
ExternalProject_Get_Property(${SAMTOOLS_PROJECT} BINARY_DIR)

SET(SAMTOOLS_INCLUDE  ${SAMTOOLS_DIR}/src/samtools_project/ CACHE INTERNAL "samtools include")
SET(SAMTOOLS_LIBRARIES
	${SAMTOOLS_DIR}/src/SAMTOOLS_PROJECT/win32/libcurses.a
	${SAMTOOLS_DIR}/src/SAMTOOLS_PROJECT/win32/libz.a
	${SAMTOOLS_DIR}/src/SAMTOOLS_PROJECT/libbam.a
	${SAMTOOLS_DIR}/src/SAMTOOLS_PROJECT/libst.a
	${HTSLIB_LIBRARY}
	CACHE INTERNAL "samtools Libraries")