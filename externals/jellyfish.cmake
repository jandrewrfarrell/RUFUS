SET_PROPERTY(DIRECTORY PROPERTY "EP_BASE" ${ep_base})

SET(JELLYFISH_PROJECT jellyfish_project CACHE INTERNAL "jellyfish project name")
SET(JELLYFISH_DIR ${CMAKE_BINARY_DIR}/externals/jellyfish CACHE INTERNAL "jellyfish project directory")
SET(JELLYFISH_LIB)



ExternalProject_Add(${JELLYFISH_PROJECT}
	GIT_REPOSITORY https://github.com/gmarcais/Jellyfish.git
	GIT_TAG master
	CONFIGURE_COMMAND "pwd; ./configure --prefix=/uufs/chpc.utah.edu/common/home/u0401321/RUFUS/bin/externals/jellyfish/src/jellyfish_project-build/"
	BUILD_IN_SOURCE 1
	#CONFIGURE_COMMAND ./configure
	BUILD_COMMAND make
	INSTALL_COMMAND make install
	UPDATE_COMMAND ""
	PREFIX ${JELLYFISH_DIR}
)

ExternalProject_Get_Property(${JELLYFISH_PROJECT} INSTALL_DIR)
ExternalProject_Get_Property(${JELLYFISH_PROJECT} SOURCE_DIR)
ExternalProject_Get_Property(${JELLYFISH_PROJECT} BINARY_DIR)

#SET(BEDTOOLS2_LIB ${SOURCE_DIR}/obj/bamToBed.o CACHE INTERNAL "GSSW Lib")
SET(JELLYFISH_INCLUDE ${SOURCE_DIR} CACHE INTERNAL "BEDTOOLS INCLUDE")