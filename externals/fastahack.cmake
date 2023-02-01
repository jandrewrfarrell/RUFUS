# Setting up external library for FASTAHACK

find_program(MAKE_EXECUTABLE NAMES make gmake mingw32-make REQUIRED)

SET(FASTAHACK_PROJECT fastahack_project CACHE INTERNAL "fastahack project name")
#SET(FASTAHACK_DIR ${CMAKE_BINARY_DIR}/externals/fastahack CACHE INTERNAL "fastahack project directory")
ExternalProject_Add(${FASTAHACK_PROJECT}
	GIT_REPOSITORY https://github.com/ekg/fastahack.git
	GIT_TAG master
	#UPDATE_DISCONNECTED true
	CONFIGURE_COMMAND ""
	BUILD_COMMAND ${MAKE_EXECUTABLE} -j -C <SOURCE_DIR>
	INSTALL_COMMAND ${MAKE_EXECUTABLE} -j -C <SOURCE_DIR> install prefix=${CMAKE_INSTALL_PREFIX}
	#BUILD_BYPRODUCTS ${my_LIBRARY}


)


ExternalProject_Get_Property(${FASTAHACK_PROJECT} INSTALL_DIR)
ExternalProject_Get_Property(${FASTAHACK_PROJECT} SOURCE_DIR)
ExternalProject_Get_Property(${FASTAHACK_PROJECT} BINARY_DIR)

SET(FASTAHACK_F_LIB ${SOURCE_DIR}/Fasta.o CACHE INTERNAL "FASTAHACK Library")
SET(FASTAHACK_S_LIB ${SOURCE_DIR}/split.o CACHE INTERNAL "FASTAHACK Library")
SET(FASTAHACK_INCLUDE ${SOURCE_DIR}/src CACHE INTERNAL "FASTAHACK Include")

#add_custom_command( MAIN_DEPENDENCY ${CMAKE_BINARY_DIR}/externals/fastahack_project-prefix/src/fastahack_project/Fasta.o
#	OUTPUT ${CMAKE_BINARY_DIR}/include/Fasta.o
#	COMMAND bash -c "mkdir -p ${CMAKE_BINARY_DIR}/include  && cp ${CMAKE_BINARY_DIR}/externals/fastahack_project-prefix/src/fastahack_project/Fasta.o ${CMAKE_BINARY_DIR}/include/"
#	VERBATIM	
#)

