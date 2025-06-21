CXX = g++
CXXFLAGS = -Iinclude -march=native -Wall -O3 -std=c++17

TARGETS = alignwithmatrix guidetree

SHARED_OBJ = align.o
UTILS_OBJ = utils.o
ALIGN_MAIN_OBJ = alignwithmatrix.o
GUIDE_MAIN_OBJ = guidetree.o

OBJS = $(SHARED_OBJ) $(UTILS_OBJ) $(ALIGN_MAIN_OBJ) $(GUIDE_MAIN_OBJ)

DEPS = include/align.hpp include/utils.hpp

all: $(TARGETS)

alignwithmatrix: $(ALIGN_MAIN_OBJ) $(SHARED_OBJ) $(UTILS_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

guidetree: $(GUIDE_MAIN_OBJ) $(SHARED_OBJ) $(UTILS_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: src/%.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGETS) $(OBJS)