CXX = g++
CXXFLAGS = -march=native -Wall -O3 -std=c++17

SRCS = alignwithmatrix.cpp

TARGET = alignwithmatrix

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET)

clean:
	rm -f $(TARGET)
