CXX = g++
CXXFLAGS = -Wall -O2 -std=c++11

SRCS = localalign.cpp

TARGET = localalign

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET)

clean:
	rm -f $(TARGET)
