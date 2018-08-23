if (CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    # when using the makefile generator, use the special variable $(MAKE) to invoke make
    # this enables the jobserver to work correctly
    set(MAKE_COMMAND "$(MAKE)")
else()
	# invoke make explicitly
	# in this case, we assume the parent build system is running in parallel already so no -j flag is added
	find_program(MAKE_COMMAND NAMES make gmake)
endif()

if (INSTALL_DEPENDENCIES)
    set(SAMTOOLS_INSTALL ${MAKE_COMMAND} install prefix=${CMAKE_INSTALL_PREFIX})
else()
	set(SAMTOOLS_INSTALL "")
endif()

# build samtools
SET(SAMTOOLS_PROJECT samtools_project CACHE INTERNAL "samtools project name")
SET(SAMTOOLS_DIR ${CMAKE_BINARY_DIR}/externals/samtools CACHE INTERNAL "samtools project directory")
ExternalProject_Add(SAMTOOLS_PROJECT
    PREFIX ${SAMTOOLS_DIR}
    GIT_REPOSITORY "https://github.com/samtools/samtools.git"
    GIT_TAG a04a7dda28e6ab209d5379c307c202b83ab67197
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ${MAKE_COMMAND} HTSDIR=${HTSLIB_INCLUDE_DIR}
    INSTALL_COMMAND "${SAMTOOLS_INSTALL}"
    LOG_DOWNLOAD 0
    LOG_UPDATE 1
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_TEST 0
    LOG_INSTALL 1
    CMAKE_CACHE_ARGS
        -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}
)

SET(SAMTOOLS_INCLUDE  ${SAMTOOLS_DIR}/src/SAMTOOLS_PROJECT/ CACHE INTERNAL "samtools include")
SET(SAMTOOLS_LIBRARIES
	${SAMTOOLS_DIR}/src/SAMTOOLS_PROJECT/win32/libcurses.a
	${SAMTOOLS_DIR}/src/SAMTOOLS_PROJECT/win32/libz.a
	${SAMTOOLS_DIR}/src/SAMTOOLS_PROJECT/libbam.a
	${SAMTOOLS_DIR}/src/SAMTOOLS_PROJECT/libst.a
	${HTSLIB_LIBRARY}
	CACHE INTERNAL "samtools Libraries")