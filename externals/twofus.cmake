SET_PROPERTY(DIRECTORY PROPERTY "EP_BASE" ${ep_base})

SET(TWOFUS_PROJECT twofus_project CACHE INTERNAL "twofus project name")
SET(TWOFUS_DIR ${CMAKE_BINARY_DIR}/externals/twofus CACHE INTERNAL "twofus project directory")
SET(TWOFUS_LIB)

ExternalProject_Add(${TWOFUS_PROJECT}
        GIT_REPOSITORY https://github.com/WilliamRichards2017/TwoFus.git
        GIT_TAG master
        INSTALL_COMMAND ""
        UPDATE_COMMAND ""
        PREFIX ${TWOFUS_DIR}
    CMAKE_CACHE_ARGS
        -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}
)

ExternalProject_Get_Property(${TWOFUS_PROJECT} INSTALL_DIR)
ExternalProject_Get_Property(${TWOFUS_PROJECT} SOURCE_DIR)
ExternalProject_Get_Property(${TWOFUS_PROJECT} BINARY_DIR)

SET(TWOFUS_LIB ${SOURCE_DIR}/bin/src/libtwofus_core.a)
SET(TWOFUS_INCLUDE ${SOURCE_DIR} CACHE INTERNAL "RUFALU Include")