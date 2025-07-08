CXX = g++
CXXFLAGS = -std=c++17 -march=native -O3 -Wall -Wextra
LDFLAGS = 

OBJS = graph.o

all: graph_degree graph_clique

graph.o: graph.cpp graph.hpp
	$(CXX) $(CXXFLAGS) -c graph.cpp -o graph.o

graph_degree: graph_degree.cpp graph.o
	$(CXX) $(CXXFLAGS) graph_degree.cpp graph.o -o graph_degree $(LDFLAGS)

graph_clique: graph_clique.cpp graph.o
	$(CXX) $(CXXFLAGS) graph_clique.cpp graph.o -o graph_clique $(LDFLAGS)

clean:
	rm -f graph_degree graph_clique *.o degree_distribution.txt