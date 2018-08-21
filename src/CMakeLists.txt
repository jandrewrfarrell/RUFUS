
INCLUDE_DIRECTORIES(
	#${SAMTOOLS_INCLUDE}
)

set(CMAKE_BUILD_TYPE Debug) # uncomment this when you need a debug build

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # clang Doesnt use pthread
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")



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
	RUFUS.interpret.cpp
	RUFUS.1kg.filter.cpp
)

SET(CORE_LIB rufus_core CACH STRING)

TARGET_LINK_LIBRARIES(${CORE_LIB})
