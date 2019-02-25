SET_PROPERTY(DIRECTORY PROPERTY "EP_BASE" ${ep_base})

SET(JELLYFISH_PROJECT jellyfish_project CACHE INTERNAL "jellyfish project name")
SET(JELLYFISH_DIR ${CMAKE_BINARY_DIR}/externals/jellyfish CACHE INTERNAL "jellyfish project directory")
SET(JELLYFISH_LIB)



ExternalProject_Add(${JELLYFISH_PROJECT}
        URL  https://github.com/gmarcais/Jellyfish/releases/download/v2.2.5/jellyfish-2.2.5.tar.gz

        CONFIGURE_COMMAND ${PROJECT_SOURCE_DIR}/bin/externals/jellyfish/src/jellyfish_project/configure --prefix=${PROJECT_SOURCE_DIR}/bin/externals/jellyfish/src/jellyfish_project/
        BUILD_IN_SOURCE 1
        BUILD_COMMAND make
        INSTALL_COMMAND make install
        UPDATE_COMMAND ""
        PREFIX ${JELLYFISH_DIR}
)

ExternalProject_Get_Property(${JELLYFISH_PROJECT} INSTALL_DIR)
ExternalProject_Get_Property(${JELLYFISH_PROJECT} SOURCE_DIR)
ExternalProject_Get_Property(${JELLYFISH_PROJECT} BINARY_DIR)

SET(JELLYFISH_INCLUDE ${SOURCE_DIR} CACHE INTERNAL "JELLYFISH INCLUDE")

