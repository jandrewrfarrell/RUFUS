SET_PROPERTY(DIRECTORY PROPERTY "EP_BASE" ${ep_base})

SET(BWA_PROJECT bwa_project CACHE INTERNAL "bwa project name")
SET(BWA_DIR ${CMAKE_BINARY_DIR}/externals/bwa CACHE INTERNAL "bwa project directory")
SET(BWA_LIB)

ExternalProject_Add( ${BWA_PROJECT}
	GIT_REPOSITORY https://github.com/williamrichards2017/bwa.git
	GIT_TAG mattOnly
	CONFIGURE_COMMAND ""
	BUILD_COMMAND "make"
	INSTALL_COMMAND ""
	UPDATE_COMMAND ""
	BUILD_IN_SOURCE 1
	PREFIX ${BWA_DIR}
)

ExternalProject_Get_Property(${BWA_PROJECT} INSTALL_DIR)
ExternalProject_Get_Property(${BWA_PROJECT} SOURCE_DIR)
ExternalProject_Get_Property(${BWA_PROJECT} BINARY_DIR)

SET(BWA_LIB ${BINARY_DIR}/bwa.o CACHE INTERNAL "BWA Lib")
SET(BWA_INCLUDE ${SOURCE_DIR} CACHE INTERNAL "BWA INCLUDE")
