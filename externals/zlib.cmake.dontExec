SET_PROPERTY(DIRECTORY PROPERTY "EP_BASE" ${ep_base})

# Setting up external library for ZLIB
SET(ZLIB_PROJECT zlib_project CACHE INTERNAL "zlib project name")
SET(ZLIB_DIR ${CMAKE_BINARY_DIR}/externals/zlib CACHE INTERNAL "zlib project directory")
ExternalProject_Add(${ZLIB_PROJECT}
	GIT_REPOSITORY https://github.com/madler/zlib.git
	GIT_TAG 50893291621658f355bc5b4d450a8d06a563053d #lock in the commit id so we don't this doesn't break in the future
    INSTALL_COMMAND "make"
    UPDATE_COMMAND ""
    BUILD_COMMAND "./configure"
    CONFIGURE_COMMAND ""
    BUILD_IN_SOURCE 1
    PREFIX ${ZLIB_DIR}
    CMAKE_CACHE_ARGS
        -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}
)

ExternalProject_Get_Property(${ZLIB_PROJECT} INSTALL_DIR)
ExternalProject_Get_Property(${ZLIB_PROJECT} SOURCE_DIR)
ExternalProject_Get_Property(${ZLIB_PROJECT} BINARY_DIR)

SET(ZLIB_LIBRARY ${BINARY_DIR}/libz.a CACHE INTERNAL "ZLIB Lib")
SET(ZLIB_LIBRARY_PATH ${BINARY_DIR} CACHE INTERNAL "ZLIB Lib Path")
SET(ZLIB_INCLUDE ${SOURCE_DIR} ${BINARY_DIR} CACHE INTERNAL "ZLIB Include")

