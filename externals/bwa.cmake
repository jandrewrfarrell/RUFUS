SET_PROPERTY(DIRECTORY PROPERTY "EP_BASE" ${ep_base})

SET(BWA_PROJECT bwa_project CACHE INTERNAL "bwa project name")
SET(BWA_DIR ${CMAKE_BINARY_DIR}/externals/bwa CACHE INTERNAL "bwa project directory")

ExternalProject_Add( ${BWA_PROJECT}
	GIT_REPOSITORY https://github.com/lh3/bwa
	GIT_TAG master
	CONFIGURE_COMMAND ""
	BUILD_COMMAND "make"
	INSTALL_COMMAND ""
	UPDATE_COMMAND ""
	PREFIX ${BWA_DIR}
)

ExternalProject_Get_Property(${BWA_PROJECT} INSTALL_DIR)
ExternalProject_Get_Property(${BWA_PROJECT} SOURCE_DIR)
ExternalProject_Get_Property(${BWA_PROJECT} BINARY_DIR)

SET