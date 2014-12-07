#OBJ_DIR = ./

BIN_SOURCES = src/ModelDist2.cpp \
		src/Overlap19.cpp \
		src/OverlapRegion2.cpp \
		src/RUFUSv5.Filter.cpp \
		src/RUFUSv6.BuildHash.cpp \
		src/ReplaceQwithDinFASTQD.cpp

#BINS = $(BIN_SOURCES:.cpp=)
BINS = $(addprefix bin/,$(notdir $(BIN_SOURCES:.cpp=)))

all: $(BINS)

CXX = g++
CXXFLAGS = -std=gnu++0x -fopenmp

bin:
	mkdir -p ./bin

$(BINS): $(BIN_SOURCES) bin
	$(CXX) src/$(notdir $@).cpp -o $@ $(CXXFLAGS)

clean:
	rm -f $(BINS)

.PHONY: clean all
