#include "graph.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>

// Calculate degree distribution
std::vector<std::pair<int, int>>
calculate_degree_distribution(const Graph &graph) {
  std::unordered_map<int, int> degree_count;

  for (const auto &neighbors : graph.adjacency_list) {
    int degree = neighbors.size();
    degree_count[degree]++;
  }

  // Convert to vector and sort by degree in descending order
  std::vector<std::pair<int, int>> degree_dist;
  for (const auto &[degree, count] : degree_count) {
    degree_dist.push_back({degree, count});
  }

  std::sort(degree_dist.begin(), degree_dist.end(),
            [](const std::pair<int, int> &a, const std::pair<int, int> &b) {
              return a.first > b.first;
            });

  return degree_dist;
}

// Output degree distribution to console
void print_degree_distribution(
    const std::vector<std::pair<int, int>> &degree_dist) {
  // std::cout << "Degree distribution (degree: count):" << std::endl;

  for (size_t i = 0; i < degree_dist.size(); i++) {
    std::cout << degree_dist[i].first << ": " << degree_dist[i].second;
    if (i < degree_dist.size() - 1) {
      std::cout << ", ";
    }
  }
  std::cout << std::endl;
}

// Save degree distribution to file for plotting
void save_degree_distribution(
    const std::vector<std::pair<int, int>> &degree_dist,
    const std::string &filename) {
  std::ofstream out_file(filename);
  if (!out_file.is_open()) {
    throw std::runtime_error("Could not create output file: " + filename);
  }

  out_file << "degree\tcount" << std::endl;
  for (const auto &[degree, count] : degree_dist) {
    out_file << degree << "\t" << count << std::endl;
  }
  out_file.close();
}

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <edge_file_path>" << std::endl;
    return 1;
  }

  try {
    std::string filename = argv[1];

    // std::cout << "Reading graph from " << filename << "..." << std::endl;
    Graph graph(filename);

    graph.clean_adjacency_lists();

    std::vector<std::pair<int, int>> degree_dist =
        calculate_degree_distribution(graph);

    print_degree_distribution(degree_dist);

    std::string output_file = "degree_distribution.txt";
    save_degree_distribution(degree_dist, output_file);

    // std::cout << "\nTotal nodes: " << graph.size() << std::endl;
    // std::cout << "Degree distribution saved to " << output_file << std::endl;

  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}