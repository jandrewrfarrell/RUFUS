set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/)

INCLUDE_DIRECTORIES(
	${FASTAHACK_INCLUDE}
)

set(CMAKE_BUILD_TYPE Debug) # uncomment this when you need a debug build

if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # clang Doesnt use pthread
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -lstdc++")

ADD_DEFINITIONS(-DBOOST_FALLTHROUGH)

SET(RUFUS_UTIL rufus_util)

add_library(${RUFUS_UTIL} STATIC
	Util.cpp
)

add_library(rufus_core STATIC
	Util.cpp
	AnnotateOverlap.cpp
	ConvertFASTqD.to.FASTQ.cpp
	ModelDist.cpp
	Overlap.cpp
	OverlapRegion.cpp
	OverlapSam.cpp
	ReplaceQwithDinFASTQD.cpp
	RUFUS.Filter.cpp
	RUFUS.Build.cpp
	RUFUS.1kg.filter.cpp
)


add_executable(AnnotateOverlap
	AnnotateOverlap.cpp	
)
target_link_libraries(AnnotateOverlap
	${RUFUS_UTIL}
)

add_executable(OverlapSam
	OverlapSam.cpp
)
target_link_libraries(OverlapSam
	${RUFUS_UTIL}
)

add_executable(ConvertFASTqD.to.FASTQ
	ConvertFASTqD.to.FASTQ.cpp
)
target_link_libraries(ConvertFASTqD.to.FASTQ
	${RUFUS_UTIL}
)

add_executable(ModelDist
	ModelDist.cpp
)
target_link_libraries(ModelDist
	${RUFUS_UTIL}
)

add_executable(Overlap
	Overlap.cpp
)
target_link_libraries(Overlap
	${RUFUS_UTIL}
)

add_executable(OverlapRegion
	OverlapRegion.cpp
)
target_link_libraries(OverlapRegion
	${RUFUS_UTIL}
)

add_executable(RUFUS.Build
	RUFUS.Build.cpp
)
target_link_libraries(RUFUS.Build
	${RUFUS_UTIL}
)

add_executable(RUFUS.Filter
	RUFUS.Filter.cpp
)
target_link_libraries(RUFUS.Filter
	${RUFUS_UTIL}
)

add_executable(RUFUS.interpret
	RUFUS.interpret.cpp
)
target_link_libraries(RUFUS.interpret
	${RUFUS_UTIL}
	${FASTAHACK_LIB}
)

add_executable(ReplaceQwithDinFASTQD
	ReplaceQwithDinFASTQD.cpp	
)
target_link_libraries(ReplaceQwithDinFASTQD
	${RUFUS_UTIL}
)


SET(CORE_LIB rufus_core CACH STRING)

TARGET_LINK_LIBRARIES(${CORE_LIB})
