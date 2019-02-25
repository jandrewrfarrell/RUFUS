SET_PROPERTY(DIRECTORY PROPERTY "EP_BASE" ${ep_base})

SET(MODIFIED_JELLYFISH_PROJECT modified_jellyfish_project CACHE INTERNAL "modifiedJellyfish project name")
SET(MODIFIED_JELLYFISH_DIR ${CMAKE_BINARY_DIR}/externals/modified_jellyfish CACHE INTERNAL "modifiedJellyfish project directory")
SET(MODIFIED_JELLYFISH_LIB)



ExternalProject_Add(${MODIFIED_JELLYFISH_PROJECT}
	URL https://github.com/WilliamRichards2017/modifiedJellyfish/blob/master/modifiedJellyfish.tar.gz

        CONFIGURE_COMMAND ${PROJECT_SOURCE_DIR}/bin/externals/modified_jellyfish/src/modified_jellyfish_project/configure --prefix=${PROJECT_SOURCE_DIR}/bin/externals/modified_jellyfish/src/modified_jellyfish_project/
        BUILD_IN_SOURCE 1
        BUILD_COMMAND make
        INSTALL_COMMAND make install
        UPDATE_COMMAND ""
        PREFIX ${MODIFIED_JELLYFISH_DIR}
)

ExternalProject_Get_Property(${MODIFIED_JELLYFISH_PROJECT} INSTALL_DIR)
ExternalProject_Get_Property(${MODIFIED_JELLYFISH_PROJECT} SOURCE_DIR)
ExternalProject_Get_Property(${MODIFIED_JELLYFISH_PROJECT} BINARY_DIR)

SET(MODIFIED_JELLYFISH_INCLUDE ${SOURCE_DIR} CACHE INTERNAL "MODIFIED JELLYFISH INCLUDE")