
SET(HTSLIB_INSTALL make install prefix=${PROJECT_SOURCE_DIR}/bin/externals/htslib/src/htslib_project})



SET(HTSLIB_PROJECT htslib_project CACHE INTERNAL "htslib project name")
SET(HTSLIB_DIR ${CMAKE_BINARY_DIR}/externals/htslib CACHE INTERNAL "htslib project directory")


ExternalProject_Add(${HTSLIB_PROJECT}
    GIT_REPOSITORY https://github.com/dillonl/htslib.git
    GIT_TAG master
    BUILD_COMMAND make lib-static
    INSTALL_COMMAND ${HTSLIB_INSTALL}
    PREFIX ${HTSLIB_DIR}
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    CMAKE_CACHE_ARGS
        -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}

)


#ExternalProject_Get_Property(${HTSLIB_PROJECT} INSTALL_DIR)
#ExternalProject_Get_Property(${HTSLIB_PROJECT} SOURCE_DIR)
#ExternalProject_Get_Property(${HTSLIB_PROJECT} BINARY_DIR)

SET(HTSLIB_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/bin/externals/htslib/src/htslib_project CACHE INTERNAL "htslib include")
SET(HTSLIB_LIBRARY ${PROJECT_SOURCE_DIR}/bin/externals/htslib/src/htslib_project/libhts.a CACHE INTERNAL "htslib Library")